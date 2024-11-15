-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

local io = terralib.includec("stdio.h")
local err = require("assert")
local base = require("base")
local tmath = require("mathfuns")
local concept = require("concept")
local array = require("arraybase")
local vecbase = require("vector")
local veccont = require("vector_contiguous")
local vecblas = require("vector_blas")
local range = require("range")

local size_t = uint64

--global flag to perform boundscheck
__boundscheck__ = true

--there is a bug on macos that leads to undefined behavior for
--simd vectors of size < 64 bytes. temporary fix is to always
--allocate a buffer that's equal or greater than 64 bytes.
local simd_fix_for_macos = function(T, N)
    local nbytes = sizeof(T) * N
    if nbytes < 64 then
        N = 64 / sizeof(T)
    end
    return N
end

local SArrayRawType = function(T, Size, options)

    --check input
    assert(terralib.types.istype(T), "ArgumentError: first argument is not a valid terra type.")
    local Size = terralib.newlist(Size)
    assert(terralib.israwlist(Size) and #Size > 0, "ArgumentError: second argument should be a list denoting the size in each dimension.")
    local Length = 1 --length of array
    for i,v in ipairs(Size) do
        assert(type(v) == "number" and v % 1 == 0 and v > 0, 
            "Expected second to last argument to be positive integers.")
        Length = Length * v
    end

    -- dimension of array
    local Dimension = #Size
    --permutation denoting order of leading dimensions. default is: {D, D-1, ... , 1}
    local Perm = options and options.perm and terralib.newlist(options.perm) or array.defaultperm(Dimension)
    array.checkperm(Perm)
    --size of leading dimension
    local SizeL = Size[Perm[1]]
    
    --generate static array struct
    local Array
    if T:isprimitive() then
        local N = simd_fix_for_macos(T, Length)
        local SIMD = vector(T, N)
        local M = sizeof(SIMD) / sizeof(T)
        Array = struct{
            union {
                data: T[M]
                simd: SIMD
            }
        }
    else
        Array = struct{
            data: T[Length]
        }
    end

    --global type traits
    Array.eltype = T
    Array.ndims = Dimension
    Array.length = Length
    Array.size = Size
    Array.perm = Perm
    Array.ldim = SizeL

    return Array
end

local SArrayStackBase = function(Array)

    local T = Array.eltype
    local N = Array.ndims
    local __size = terralib.constant(terralib.new(size_t[N], Array.size))
    
    terra Array:length() : size_t
        return [Array.length]
    end

    terra Array:getdataptr() : &T
        return &self.data[0]
    end

    --element size as a terra method
    Array.methods.size = terralib.overloadedfunction("size", {
        terra(self : &Array, i : size_t)
            return __size[i]
        end,
        terra(self : &Array)
            return __size
        end
    })

    --get lowlevel base functionality for nd arrays
    local arraybase = array.ArrayBase(Array)

    --dimension permutation as a terra method
    Array.methods.perm = arraybase.getperm

    --given multi-indices {i_1, i_2, ..., i_D}
    --compute the linear index as follows (assuming perm={D, D-1, ... , 1}):
    --l = i_d + 
    --    Size[D] * i_{d-1} + 
    --    Size[D] * Size[D-1] * i_{d-2} +
    --    .
    --    .
    --    Size[D] * Size[D-1] * ... * Size[2] * i_{1}
    local Indices = arraybase.Indices

    local getlinearindex
    if N == 1 then
        terra getlinearindex(self : &Array, index : size_t) : size_t
            arraybase.boundscheck(self, index)
            return index
        end
    else
        terra getlinearindex(self : &Array, [Indices]) : size_t
            arraybase.boundscheck(self, Indices)
            var lindex = [ Indices[ Array.perm[1] ] ]
            escape
                local cumprod = Array.size[ Array.perm[1] ]
                for k = 1, N-1 do
                    emit quote
                        lindex = lindex + cumprod * [ Indices[ Array.perm[k+1] ] ]
                    end
                    cumprod = cumprod * Array.size[ Array.perm[k+1] ]
                end
            end
            return lindex
        end
    end
    getlinearindex:setinlined(true)

    local get = terra(self : &Array, index : size_t)
        arraybase.boundscheck(self, index)
        return self.data[index]
    end
    
    local set = terra(self : &Array, index : size_t, x : T)
        arraybase.boundscheck(self, index)
        self.data[index] = x
    end
    
    if N == 1 then
        Array.methods.get = get
        Array.methods.set = set
    else
        Array.methods.get = terralib.overloadedfunction("get", 
        {
            get,
            terra(self : &Array, [Indices])
                return self.data[getlinearindex(self, [Indices])]
            end
        })
        
        Array.methods.set = terralib.overloadedfunction("set",
        {
            set, 
            terra(self : &Array, [Indices], x : T)
                self.data[getlinearindex(self, [Indices])] = x
            end
        })
    end

    Array.metamethods.__apply = macro(function(self, ...)
        local indices = terralib.newlist{...}
        if #indices == 1 then
            if indices[1].tree.type.convertible == "tuple" then
                if self:gettype():ispointer() then
                    return `self.data[getlinearindex(self, unpacktuple([indices]))]
                else
                    return `self.data[getlinearindex(&self, unpacktuple([indices]))]
                end
            else
                local index = indices[1] 
                return `self.data[ [index] ]
            end
        else   
            if self:gettype():ispointer() then
                return `self.data[getlinearindex(self, [indices])]
            else
                return `self.data[getlinearindex(&self, [indices])]
            end
        end
    end)

    Array.staticmethods.new = terra()
        return Array{}
    end

    local checkarraysize = function (arraysize)
        for k = 1, N do
            assert(arraysize[k] == Array.size[k], "ArgumentError: sizes in dimension " .. tostring(k) .. " is not consistent with array dimensions.")
        end
    end

    Array.staticmethods.from = macro(function(args)
        local arrayentries = terralib.newlist(args:asvalue())
        local arraysize = arraybase.processarrayinput(arrayentries)
        checkarraysize(arraysize)
        return quote
            var A : Array
            escape
                for mi in array.productiter(unpack(arraysize)) do
                    emit quote 
                        A:set([ mi ], [ array.getarrayentry(args, mi) ]) 
                    end
                end
            end
        in
            A
        end
    end)
    
end

local SArrayVectorBase = function(Array)

    local T = Array.eltype

    vecbase.VectorBase(Array) --add fall-back routines

    Array.staticmethods.all = terra(value : T)
        var A : Array
        for i = 0, A:length() do
            A:set(i, value)
        end
        return A
    end

    if concept.Number(T) then

        Array.staticmethods.zeros = terra()
            return Array.all(T(0))
        end

        Array.staticmethods.ones = terra()
            return Array.all(T(1))
        end

    end

    --specializations using simd operations
    if T:isprimitive() then

        Array.staticmethods.all = terra(value : T)
            var A = Array.new()
            A.simd = value
            return A
        end

        terra Array:fill(v : T)
            self.simd = v
        end

        terra Array:copy(other : &Array)
            self.simd = other.simd
        end

        terra Array:scal(a : T)
            self.simd = a * self.simd
        end

        terra Array:axpy(a : T, x : &Array)
            self.simd = self.simd + a * x.simd
        end

    end

end

local SArrayIteratorBase = function(Array)

    local T = Array.eltype
    local N = Array.ndims
    local Unitrange = range.Unitrange(int)

    local __uranges = Array.size:map(function(s) return terralib.constant( terralib.new(Unitrange, {0, s}) ) end)

    --get lowlevel base functionality for nd arrays
    local arraybase = array.ArrayBase(Array)

    --return linear indices product range
    terra Array:linear_indices()
        return Unitrange{0, [Array.length]}
    end

    --return linear indices product range
    terra Array:unitrange(i : size_t)
        return Unitrange{0, self:size(i)}
    end

    --return cartesian indices product range
    terra Array:cartesian_indices()
        return range.product([__uranges], {perm = {[Array.perm]}})
    end

    terra Array:rowmajor_cartesian_indixes()
        return range.product([__uranges], {perm = {[array.defaultperm(N)]}})
    end

    --printing array contents
    Array.methods.print = arraybase.print

    --standard iterator is added in VectorBase
    vecbase.IteratorBase(Array) --add fall-back routines

end

local StaticArray = function(T, Size, options)
    
    --generate the raw type
    local Array = SArrayRawType(T, Size, options)
    
    --print typename
    function Array.metamethods.__typename(self)
        local sizes = "{"
        local perm = "{"
        for i = 1, Array.ndims-1 do
            sizes = sizes .. tostring(Array.size[i]) .. ","
            perm = perm .. tostring(Array.perm[i]) .. ","
        end
        sizes = sizes .. tostring(Array.size[Array.ndims]) .. "}"
        perm = perm .. tostring(Array.perm[Array.ndims]) .. "}"
        return "StaticArray(" .. tostring(T) ..", " .. sizes .. ", perm = " .. perm .. ")"
    end
    
    --add base functionality
    base.AbstractBase(Array)

    --implement interfaces
    SArrayStackBase(Array)
    SArrayVectorBase(Array)
    SArrayIteratorBase(Array)

    return Array
end

local StaticVector = terralib.memoize(function(T, N)
    
    --generate the raw type
    local SVector = SArrayRawType(T, {N})
    
    function SVector.metamethods.__typename(self)
        return ("StaticVector(%s, %d)"):format(tostring(T), N)
    end
    
    --add base functionality
    base.AbstractBase(SVector)

    --implement interfaces
    SArrayStackBase(SVector)
    SArrayVectorBase(SVector)
    SArrayIteratorBase(SVector)

    veccont.VectorContiguous:addimplementations{SVector}

    if concept.BLASNumber(T) then
        terra SVector:getblasinfo()
            return self:length(), self:getdataptr(), 1
        end
        vecblas.VectorBLAS:addimplementations{SVector}
    end

    return SVector
end)

local StaticMatrix = terralib.memoize(function(T, Size, options)
    
    local SMatrix = SArrayRawType(T, Size, options)

    --check that a matrix-type was generated
    assert(SMatrix.ndims == 2, "ArgumentError: second argument should be a table with matrix dimensions.")

    function SVector.metamethods.__typename(self)
        return ("StaticMatrix(%s, {%d, %d})"):format(tostring(T), Size{1}, Size{2})
    end

    --add base functionality
    base.AbstractBase(SMatrix)

    --implement interfaces
    SArrayStackBase(SMatrix)
    SArrayVectorBase(SMatrix)
    SArrayIteratorBase(SMatrix)

    if concept.BLASNumber(T) then
        terra SMatrix:getblasdenseinfo()
            return [ SMatrix.size[1] ], [ SMatrix.size[2] ], self:getdataptr(), [ SMatrix.ldim ]
        end
        local matblas = require("matrix_blas_dense")
        matblas.BLASDenseMatrixBase(SMatrix)
    end

    return SMatrix
end)


return {
    StaticArray = StaticArray,
    SArrayRawType = SArrayRawType,
    SArrayStackBase = SArrayStackBase,
    SArrayVectorBase = SArrayVectorBase,
    SArrayIteratorBase = SArrayIteratorBase,
    StaticVector = StaticVector,
    StaticMatrix = StaticMatrix
}