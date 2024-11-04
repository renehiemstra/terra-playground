-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

import "terraform"
local io = terralib.includec("stdio.h")
local err = require("assert")
local base = require("base")
local vecbase = require("vector")
local veccont = require("vector_contiguous")
local vecblas = require("vector_blas")
local concept = require("concept")
local rn = require("range")

local size_t = uint64

local ntuple = function(T,N)
    assert(terralib.types.istype(T), "ArgumentError: Not a valid terra type.")
    local types = terralib.newlist()
    for k = 1, N do
        types:insert(T)
    end
    return tuple(unpack(types))
end

--global flag to perform boundscheck
__boundscheck__ = true

local StaticArray = terralib.memoize(function(T, ...)
    --check input
    local Size = terralib.newlist{...}
    local M = 1
    for i,v in ipairs(Size) do
        assert(type(v) == "number" and v % 1 == 0 and v > 0, 
            "Expected second to last argument to be positive integers.")
        M = M * v
    end
    local D = #Size
    local N = Size[D] --leading dimension D
    local SIMD = vector(T, N) --perform simd operations along dimensions D

    --define struct for static array, stored as a vector with
    --set and get methods for a tensor of type defined by its Size table
    local sarray = struct{
        union {
            data: T[M]
            simd: SIMD
        }
    }
    --print typename
    function sarray.metamethods.__typename(self)
        local sizes = "{"
        for i = 1, D-1 do
            sizes = sizes .. tostring(Size[i]) .. ","
        end
        sizes = sizes .. tostring(Size[D]) .. "}"
        return "StaticArray(" .. tostring(T) .."," .. sizes .. ")"
    end
    --add base functionality
    base.AbstractBase(sarray)

    --global type traits
    sarray.eltype = T
    --local type traits
    local __size = terralib.constant(terralib.new(size_t[D], Size))

    --element size as a static lua method
    function sarray.size(k)
        return k and Size[k] or Size
    end

    --element size as a terra method
    sarray.methods.size = terralib.overloadedfunction("size", {
        terra(self : &sarray, i : size_t)
            return __size[i]
        end,
        terra(self : &sarray)
            return __size
        end
    })

    terra sarray:length()
        return M
    end

    --given multi-indices {i_1, i_2, ..., i_D}
    --compute the linear index as follows:
    --l = i_d + 
    --    Size[D] * i_{d-1} + 
    --    Size[D] * Size[D-1] * i_{d-2} +
    --    .
    --    .
    --    Size[D] * Size[D-1] * ... * Size[2] * i_{1}
    local Indices = terralib.newlist{}
    for i = 1, D do
        Indices:insert(symbol(size_t))
    end

    local terra __boundscheck([Indices])
        escape
            for d = 1, D do
                local index, size = Indices[d], Size[d]
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
            return `__boundscheck([indices])
        end
    end)

    local terra getlinearindex([Indices])
        boundscheck(Indices)
        var lindex = [ Indices[D] ]
        escape
            local cumprod = Size[D]
            for d = D, 2, -1 do
                emit quote
                    lindex = lindex + cumprod * [ Indices[d-1] ]
                end
                cumprod = cumprod * Size[d-1]
            end
        end
        return lindex
    end
    getlinearindex:setinlined(true)

    terra sarray:get([Indices])
        var lindex = getlinearindex([Indices])
        return self.data[lindex]
    end
    
    terra sarray:set([Indices], x : T)
        var lindex = getlinearindex([Indices])
        self.data[lindex] = x
    end

    sarray.metamethods.__apply = macro(function(self, ...)
        local indices = terralib.newlist{...}
        return quote
            var lindex = getlinearindex([indices])
        in
            self.data[lindex]
        end
    end)

    return sarray
end)


import "terratest/terratest"

testenv "Arbitrary dimension arrays" do

    local linrange = rn.Unitrange(int)
    local Array = StaticArray(int, 2, 3, 4)

    testset "3D arrays" do
        terracode
            var A : Array
            for count, indices in rn.enumerate(rn.product(linrange{0,2}, linrange{0,3}, linrange{0,4})) do
                A:set(unpacktuple(indices), count)
            end
        end
        test A:size(0)  == 2
        test A:size(1)  == 3
        test A:size(2)  == 4
        test A:length() == 24 
        for k = 0, 23 do
            test A.data[k] == k
        end
    end
end





--return {
----    StaticArray = StaticArray
--}
