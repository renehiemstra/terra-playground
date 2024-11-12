-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

local io = terralib.includec("stdio.h")
local err = require("assert")
local base = require("base")
local tmath = require("mathfuns")
local concept = require("concept")
local vecbase = require("vector")
local range = require("range")

local size_t = uint64


--global flag to perform boundscheck
__boundscheck__ = true

local function checkperm(perm)
    assert(terralib.israwlist(perm), "ArgumentError: input should be a raw list.")
    local linrange = perm:mapi(function(i,v) return i end)
    for i,v in ipairs(perm) do
        linrange[v] = nil
    end
    assert(#linrange == 0, "ArgumentError: input list is not a valid permutation.")
end

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
    local Perm = options and options.perm and terralib.newlist(options.perm) or Size:mapi(function(i,v) return Dimension+1-i end)
    checkperm(Perm)
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




local linrange = {}
linrange.__index = linrange
linrange.__metatable = linrange

function linrange.new(n)
    local t = {}
    t.n = n
    return setmetatable(t, linrange)
end

function linrange:next(i)
    local i = i+1
    if i < self.n then
        return i
    end
end

function linrange:iterate()
    return linrange.next, self, -1
end

local function next(sizes, iter, i, k)
    i[k] = iter[k]:next(i[k])
    if i[k] then
        return true
    else
        if k>1 then
            iter[k] = linrange.new(sizes[k])
            i[k] = iter[k]:next(-1)
            return next(sizes, iter, i, k-1)
        end
    end
    return false
end

local productiter = function(...)
    local sizes = terralib.newlist{...}
    local m = #sizes
    local iter = sizes:map(function(s) return linrange.new(s) end)
    local i = sizes:map(function(s) return 0 end)
    i[m] = -1
    return function()
        if next(sizes, iter, i, m) then
            return i
        end
    end
end

local SArrayStackBase = function(Array)

    local T = Array.eltype
    local __size = terralib.constant(terralib.new(size_t[Array.ndims], Array.size))
    local __perm = terralib.constant(terralib.new(size_t[Array.ndims], Array.perm))

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

    --element size as a terra method
    Array.methods.perm = terralib.overloadedfunction("perm", {
        terra(self : &Array, i : size_t)
            return __perm[i]
        end,
        terra(self : &Array)
            return __perm
        end
    })

    --given multi-indices {i_1, i_2, ..., i_D}
    --compute the linear index as follows (assuming perm={D, D-1, ... , 1}):
    --l = i_d + 
    --    Size[D] * i_{d-1} + 
    --    Size[D] * Size[D-1] * i_{d-2} +
    --    .
    --    .
    --    Size[D] * Size[D-1] * ... * Size[2] * i_{1}
    local Indices = Array.size:map(function(s) return symbol(size_t) end)

    local terra __boundscheck_linear(index : size_t)
        err.assert(index < [ Array.length ], "BoundsError: out of bounds.")
    end

    local terra __boundscheck_cartesian([Indices])
        escape
            for d = 1, Array.ndims do
                local index, size = Indices[d], Array.size[d]
                local message = "BoundsError: array dimension " .. tostring(d) .. " out of bounds."
                emit quote
                    err.assert([index] < [size], message)
                end
            end
        end
    end

    local boundscheck = macro(function(...)
        if __boundscheck__ then
            local indices = terralib.newlist{...}
            if #indices == 1 then
                return `__boundscheck_linear([indices])
            else
                return `__boundscheck_cartesian([indices])
            end
        end
    end)

    local getlinearindex
    if Array.ndims == 1 then
        terra getlinearindex(index : size_t) : size_t
            boundscheck(index)
            return index
        end
    else
        terra getlinearindex([Indices]) : size_t
            boundscheck(Indices)
            var lindex = [ Indices[ Array.perm[1] ] ]
            escape
                local cumprod = Array.size[ Array.perm[1] ]
                for k = 1, Array.ndims-1 do
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
        boundscheck(index)
        return self.data[index]
    end
    
    local set = terra(self : &Array, index : size_t, x : T)
        boundscheck(index)
        self.data[index] = x
    end
    
    if Array.ndims == 1 then
        Array.methods.get = get
        Array.methods.set = set
    else
        Array.methods.get = terralib.overloadedfunction("get", 
        {
            get,
            terra(self : &Array, [Indices])
                return self.data[getlinearindex([Indices])]
            end
        })
        
        Array.methods.set = terralib.overloadedfunction("set",
        {
            set, 
            terra(self : &Array, [Indices], x : T)
                self.data[getlinearindex([Indices])] = x
            end
        })
    end

    Array.metamethods.__apply = macro(function(self, ...)
        local indices = terralib.newlist{...}
        if #indices == 1 then
            if indices[1].tree.type.convertible == "tuple" then
                return `self.data[getlinearindex(unpacktuple([indices]))]
            else
                local index = indices[1] 
                return `self.data[ [index] ]
            end
        else   
            return `self.data[getlinearindex([indices])]
        end
    end)

    Array.staticmethods.new = terra()
        return Array{}
    end

    local function processarraydim(array, k)
        assert(#array == Array.size[k], "ArgumentError: array input size inconsistent with array dimensions.")
        if k < Array.ndims then
            for i,v in ipairs(array) do
                assert(type(v) == "table", "ArgumentError: expected array input of dimension " .. tostring(Array.ndims) ..".")
                array[i] = terralib.newlist(v)
                assert(#array[i] == Array.size[k+1], "ArgumentError: array input size inconsistent with array dimensions.")
                processarraydim(array[i], k+1)
            end
        end
    end

    local function processarrayinput(array)
        processarraydim(array, 1)
        return array
    end

    local function getarrayentry(array, mi, k)
        local k = k or 1
        if type(array) == "table" then
            return getarrayentry(array[ mi[k]+1 ], mi, k+1)
        else
            return array
        end
    end

    Array.staticmethods.from = macro(function(args)
        local arrayentries = terralib.newlist(args:asvalue())
        arrayentries = processarrayinput(arrayentries)
        return quote
            var array : Array
            escape
                for mi in productiter(unpack(Array.size)) do
                    local v = getarrayentry(arrayentries, mi)
                    emit quote array([mi]) = [ v ] end
                end
            end
        in
            array
        end
    end)
    
end

local SArrayVectorBase = function(Array)

    local T = Array.eltype

    terra Array:length()
        return [ Array.length ]
    end

    vecbase.VectorBase(Array) --add fall-back routines

    Array.staticmethods.all = terra(value : T)
        var A = Array.new()
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

    local Unitrange = range.Unitrange(int)
    local __uranges = Array.size:map(function(s) return terralib.constant( terralib.new(Unitrange, {0, s}) ) end)
    local rowmajorperm = Array.size:mapi(function(i,v) return Array.ndims+1-i end)

    --return linear indices product range
    terra Array:linear_indices()
        return Unitrange{0, [Array.length]}
    end

    --return cartesian indices product range
    terra Array:cartesian_indices()
        return range.product([__uranges], {perm = {[Array.perm]}})
    end

    terra Array:rowmajor_cartesian_indixes()
        return range.product([__uranges], {perm = {[rowmajorperm]}})
    end

    --standard iterator is added in VectorBase
    vecbase.IteratorBase(Array) --add fall-back routines

    local print_vector = function(self, name)
        return quote
            io.printf("%s = \n", [ name ])
            for v in self do
                io.printf("[%s]\n", tmath.numtostr(v))
            end
            io.printf("\n")
        end
    end

    local print_matrix = function(self, name)
        local I, J = __uranges[1], __uranges[2]
        return quote 
            io.printf("%s = \n", [ name ])
            for i in I do
                io.printf("\t[")
                for j in J do
                    var value = self(i, j)
                    io.printf("%s\t", tmath.numtostr(value))
                end
                io.printf("]\n")
            end
            io.printf("\n")
        end
    end

    local print_array = function(self, name)
        local m = Array.ndims
        local I, J = __uranges[m-1], __uranges[m]
        local K = __uranges:filteri(function(i,v) return i <= m - 2 end)
        local p = rowmajorperm:filteri(function(i,v) return i <= m - 2 end)
        local ntimes = K:mapi(function(i,v) return "%d" end)
        local slice = name .."[" .. table.concat(ntimes,",") .. ", :, :] = \n"
        return quote 
            for k in range.product([K]) do
                io.printf([ slice ], unpacktuple(k))
                for i in I do
                    io.printf("\t[")
                    for j in J do
                        var value = self(unpacktuple(k), i, j)
                        io.printf("%s\t", tmath.numtostr(value))
                    end
                    io.printf("]\n")
                end
                io.printf("\n")
            end
        end
    end

    Array.methods.show = macro(function(self)
        local name = self.tree.name
        if Array.ndims == 1 then
            return print_vector(self, name)
        elseif Array.ndims == 2 then
            return print_matrix(self, name)
        else
            return print_array(self, name)
        end
    end)
 
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


return {
    StaticArray = StaticArray,
    SArrayRawType = SArrayRawType,
    SArrayStackBase = SArrayStackBase,
    SArrayVectorBase = SArrayVectorBase,
    SArrayIteratorBase = SArrayIteratorBase
}