-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

local io = terralib.includec("stdio.h")
local err = require("assert")
local alloc = require("alloc")
local base = require("base")
local tmath = require("mathfuns")
local concept = require("concept")
local array = require("arraybase")
local vecbase = require("vector")
local matbase = require("matrix")
local veccont = require("vector_contiguous")
local vecblas = require("vector_blas")
local tup = require("tuple")
local range = require("range")

local luafun = require("fun")

local Allocator = alloc.Allocator
local size_t = uint64

--global flag to perform boundscheck
__boundscheck__ = true


local DArrayRawType = function(T, Dimension, options)

    --check input
    assert(terralib.types.istype(T), "ArgumentError: first argument is not a valid terra type.")

    --permutation denoting order of leading dimensions. default is: {D, D-1, ... , 1}
    local Perm = options and options.perm and terralib.newlist(options.perm) or array.defaultperm(Dimension)
    array.checkperm(Perm)
    
    --generate static array struct
    local S = alloc.SmartBlock(T)

    local struct Array{
        data : S
        size : size_t[Dimension]
        cumsize : size_t[Dimension] --cumulative product dimensions - is ordered according to 'perm'
    }

    --global type traits
    Array.eltype = T
    Array.ndims = Dimension
    Array.perm = Perm

    return Array
end

local DArrayStackBase = function(Array)

    local T = Array.eltype
    local N = Array.ndims
    local Sizes = tup.ntuple(size_t, N)

    terra Array:length() : size_t
        return self.cumsize[ [N-1] ]
    end

    terra Array:getdataptr() : &T
        return self.data:getdataptr()
    end

    --element size as a terra method
    Array.methods.size = terralib.overloadedfunction("size", {
        terra(self : &Array, i : size_t)
            return self.size[i]
        end,
        terra(self : &Array)
            return self.size
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
                for k = 1, N-1 do
                    emit quote
                        lindex = lindex + self.cumsize[ [k-1] ] * [ Indices[ Array.perm[k+1] ] ]
                    end
                end
            end
            return lindex
        end
    end
    getlinearindex:setinlined(true)

    local get = terra(self : &Array, index : size_t)
        arraybase.boundscheck(self, index)
        return self.data:get(index)
    end
    
    local set = terra(self : &Array, index : size_t, x : T)
        arraybase.boundscheck(self, index)
        self.data:set(index, x)
    end
    
    if N == 1 then
        Array.methods.get = get
        Array.methods.set = set
    else
        Array.methods.get = terralib.overloadedfunction("get", 
        {
            get,
            terra(self : &Array, [Indices])
                return self.data:get(getlinearindex(self, [Indices]))
            end
        })
        
        Array.methods.set = terralib.overloadedfunction("set",
        {
            set, 
            terra(self : &Array, [Indices], x : T)
                self.data:set(getlinearindex(self, [Indices]), x)
            end
        })
    end

    Array.metamethods.__apply = macro(function(self, ...)
        local indices = terralib.newlist{...}
        if #indices == 1 then
            local index = indices[1]
            if index.tree.type.convertible == "tuple" then
                if self:gettype():ispointer() then
                    return `self.data( getlinearindex(self, unpacktuple([index])) )
                else
                    return `self.data( getlinearindex(&self, unpacktuple([index])) )
                end
            else
                return `self.data( [index] )
            end
        else
            if self:gettype():ispointer() then
                return `self.data( getlinearindex(self, [indices]) )
            else
                return `self.data( getlinearindex(&self, [indices]) )
            end
        end
    end)

    local terra getcumsize(size : size_t[N])
        var cumsize : size_t[N]
        escape
            local p = Array.perm[1]
            emit quote cumsize[0] = size[ [p-1] ] end
            for k = 2, N do
                local p = Array.perm[k]
                emit quote cumsize[ [k-1] ] = cumsize[ [k-2] ] * size[ [p - 1] ] end
            end
        end
        return cumsize
    end

    --create a new dynamic array
    local new = terra(alloc: Allocator, size : tup.ntuple(size_t, N))
        var __size = [ &size_t[N] ](&size)  --we need the size as an array
        var cumsize = getcumsize(@__size)   --compute cumulative sizes
        var length = cumsize[N-1]           --length is last entry in 'cumsum'
        return Array{alloc:allocate(sizeof(T), length), @__size, cumsize}
    end

    --For N==1 we allow passing the size as an integer or as a tuple holding
    --a single integer
    if N==1 then
        Array.staticmethods.new = terralib.overloadedfunction("new", {
            new,
            terra(alloc: Allocator, size : size_t)
                return new(alloc, {size})
            end
        })
    else
        Array.staticmethods.new = new
    end

    Array.staticmethods.from = macro(function(alloc, args)
        local arrayentries = terralib.newlist(args:asvalue())
        local arraysize = arraybase.processarrayinput(arrayentries)
        return quote
            var A = Array.new([alloc], {[arraysize]})
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

    terra Array:similar()
    end

end

local DArrayVectorBase = function(Array)

    local T = Array.eltype
    local N = Array.ndims
    local Sizes = tup.ntuple(size_t, N)

    vecbase.VectorBase(Array) --add fall-back routines

    local all = function(S)
        return terra(alloc : Allocator, size : S, value : T)
            var A = Array.new(alloc, size)
            for i = 0, A:length() do
                A:set(i, value)
            end
            return A
        end
    end

    local zeros = function(S)
        return terra(alloc : Allocator, size : S)
            return Array.all(alloc, size, T(0))
        end
    end

    local ones = function(S)
        return terra(alloc : Allocator, size : S)
            return Array.all(alloc, size, T(1))
        end
    end


    if N == 1 then

        Array.staticmethods.all = terralib.overloadedfunction("all", {all(size_t), all(Sizes)})

        if concept.Number(T) then
            Array.staticmethods.zeros = terralib.overloadedfunction("zeros", {zeros(size_t), zeros(Sizes)})
            Array.staticmethods.ones = terralib.overloadedfunction("ones", {ones(size_t), ones(Sizes)})
        end

    else
        Array.staticmethods.all = all(Sizes)
        
        if concept.Number(T) then
            Array.staticmethods.zeros = zeros(Sizes)
            Array.staticmethods.ones = ones(Sizes)
        end
    end

end

local DArrayIteratorBase = function(Array)

    local T = Array.eltype
    local N = Array.ndims
    local Unitrange = range.Unitrange(int)

    --get lowlevel base functionality for nd arrays
    local arraybase = array.ArrayBase(Array)

    --return linear indices product range
    terra Array:linear_indices()
        return Unitrange{0, self:length()}
    end

    --return linear indices product range
    terra Array:unitrange(i : size_t)
        return Unitrange{0, self:size(i)}
    end

    --get unitranges in all array dimensions
    local unitranges = arraybase.getunitranges(N)

    --return cartesian indices product range
    terra Array:cartesian_indices()
        var K = unitranges(self)
        return range.product(unpacktuple(K), {perm = {[Array.perm]}})
    end

    terra Array:rowmajor_cartesian_indixes()
        var K = unitranges(self)
        return range.product(unpacktuple(K), {perm = {[array.defaultperm(N)]}})
    end

    --printing array contents
    Array.methods.print = arraybase.print

    --standard iterator is added in VectorBase
    vecbase.IteratorBase(Array) --add fall-back routines
 
end

local DynamicArray = function(T, Dimension, options)
    
    --generate the raw type
    local Array = DArrayRawType(T, Dimension, options)
    
    --print typename
    function Array.metamethods.__typename(self)
        local sizes = "{"
        local perm = "{"
        for i = 1, Array.ndims-1 do
            perm = perm .. tostring(Array.perm[i]) .. ","
        end
        perm = perm .. tostring(Array.perm[Array.ndims]) .. "}"
        return "DynamicArray(" .. tostring(T) ..", " .. tostring(Array.ndims) .. ", perm = " .. perm .. ")"
    end
    
    --add base functionality
    base.AbstractBase(Array)

    --implement interfaces
    DArrayStackBase(Array)
    DArrayVectorBase(Array)
    DArrayIteratorBase(Array)

    return Array
end

local DynamicVector = function(T)

    --generate the raw type
    local DVector = DArrayRawType(T, 1)
    
    function DVector.metamethods.__typename(self)
        return ("DynamicVector(%s)"):format(tostring(T))
    end
    
    --add base functionality
    base.AbstractBase(DVector)

    --implement interfaces
    DArrayStackBase(DVector)
    DArrayVectorBase(DVector)
    DArrayIteratorBase(DVector)

    veccont.VectorContiguous:addimplementations{DVector}

    if concept.BLASNumber(T) then
        terra DVector:getblasinfo()
            return self:length(), self:getdataptr(), 1
        end
        vecblas.VectorBLAS:addimplementations{DVector}
    end

    assert(vecbase.Vector(DVector))

    return DVector
end

local DynamicMatrix = function(T, options)
    
    local DMatrix = DArrayRawType(T, 2, options)

    --check that a matrix-type was generated
    assert(DMatrix.ndims == 2, "ArgumentError: expected an array of dimension 2.")

    function DMatrix.metamethods.__typename(self)
        local perm = "{"
        for i = 1, Array.ndims-1 do
            perm = perm .. tostring(Array.perm[i]) .. ","
        end
        perm = perm .. tostring(Array.perm[Array.ndims]) .. "}"
        return "DynamicMatrix(" .. tostring(T) ..", " .. tostring(Array.ndims) .. ", perm = " .. perm .. ")"
    end

    --add base functionality
    base.AbstractBase(DMatrix)

    --implement interfaces
    DArrayStackBase(DMatrix)
    DArrayVectorBase(DMatrix)
    DArrayIteratorBase(DMatrix)

    --add linear operator functionality
    matbase.MatrixBase(DMatrix)
    
    if concept.BLASNumber(T) then
        terra DMatrix:getblasdenseinfo()
            return self:size(0), self:size(1), self:getdataptr(), self:size(1)
        end
        local matblas = require("matrix_blas_dense")
        matblas.BLASDenseMatrixBase(DMatrix)
    end

    assert(matbase.Matrix(DMatrix))

    return DMatrix
end


return {
    DArrayRawType = DArrayRawType,
    DArrayStackBase = DArrayStackBase,
    DArrayIteratorBase = DArrayIteratorBase,
    DynamicArray = DynamicArray,
    DynamicVector = DynamicVector,
    DynamicMatrix = DynamicMatrix
}