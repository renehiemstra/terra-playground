-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

local err = require("assert")
local base = require("base")
local vecbase = require("vector")
local veccont = require("vector_contiguous")
local vecblas = require("vector_blas")
local concept = require("concept")

local StaticVector = terralib.memoize(function(T, N)
    local function create_static_vector(T, N)
        local V
        if T:isprimitive() then
            local nbytes = sizeof(T) * N
            if nbytes < 64 then
                N = 64 / sizeof(T)
            end
            local SIMD = vector(T, N)
            local M = sizeof(SIMD) / sizeof(T)
            V = struct{
                union {
                    data: T[M]
                    simd: SIMD
                }
            }
        else
           V = struct{
                data: T[N]
            }
        end
        function V.metamethods.__typename(self)
            return ("StaticVector(%s, %d)"):format(tostring(T), N)
        end
        base.AbstractBase(V)
        return V
    end

    local V = create_static_vector(T, N)
    V.eltype = T

    terra V:size(): uint64
        return N
    end

    terra V:get(i: uint64)
        err.assert(i < N)
        return self.data[i]
    end

    terra V:set(i: uint64, x: T)
        err.assert(i < N)
        self.data[i] = x
    end

    V.metamethods.__apply = macro(function(self, i)
        return quote err.assert(i < N) in self.data[i] end
    end)

    terra V:getbuffer()
        return self:size(), &self.data[0]
    end

    veccont.VectorContiguous:addimplementations{V}

    terra V:getblasinfo()
        return self:size(), &self.data, 1
    end

    vecblas.VectorBLAS:addimplementations{V}

    V.staticmethods.new = terra()
        return V {}
    end

    V.staticmethods.all = terra(value: T)
        var v = V.new()
        escape
            for i = 0, N - 1 do
                emit quote v.data[i] = value end
            end
        end
        return v
    end

    V.staticmethods.zeros = terra()
        return V.all(0)
    end

    V.staticmethods.ones = terra()
        return V.all(1)
    end

    V.staticmethods.from = macro(
        function(...)
            local args = {...}
            assert(#args == N, "Length of input list does not match static dimension")
            local vec = symbol(V)
            local set_values = terralib.newlist()
            for i, v in ipairs(args) do
                set_values:insert(quote [vec].data[i - 1] = v end)
            end
            return quote
                var [vec] = V.new()
                [set_values]
            in
                [vec]
            end
        end
    )

    vecbase.VectorBase(V)

    if T:isprimitive() then
        V.templates.fill[{&V.Self, concept.Number} -> {}] = function(Self, S)
            local terra fill(self: Self, a: S)
                self.simd = a
            end
            return fill
        end

        V.templates.copy[{&V.Self, &V.Self} -> {}] = function(V1, V2)
            local terra copy(self: V1, other: V2)
                self.simd = other.simd
            end
            return copy
        end

        V.templates.scal[{&V.Self, concept.Number} -> {}] = function(Self, S)
            local terra scal(self: Self, a: S)
                self.simd = a * self.simd
            end
            return scal
        end

        V.templates.axpy[{&V.Self, concept.Number, &V.Self} -> {}] =
            function(V1, S, V2)
                local terra axpy(y: V1, a: S, x: V2)
                    y.simd = y.simd + a * x.simd
                end
                return axpy
            end
        -- dot impletation doesn't profit from a vectorized implementation
        -- as the operation works vertically and thus requires synchronization
        -- of the vector registers.
    end

    return V
end)

return {
    StaticVector = StaticVector
}
