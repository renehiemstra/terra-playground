-- SPDX-FileCopyrightText: 2025 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2025 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT
local alloc = require("alloc")
local complex = require("complex")
local contraction = require("contraction")
local darray = require("darray")
local nfloat = require("nfloat")
local tmath = require("tmath")

import "terratest/terratest"

local complex128 = complex.complex(double)
local nfloat256 = nfloat.FixedFloat(256)

local types = {
    [float] = `1e-6f,
    [double] = `1e-15,
    [complex128] = `1e-14,
    [nfloat256] = `1e-60,
}

for T, tol in pairs(types) do
    testenv(T) "Symmetric bilinear contraction" do
        terracode
            var A: alloc.DefaultAllocator()
            var nx = 1
            var nt = 51
            var nv = 25
            var q = [darray.DynamicArray(T, 3)].new(&A, {nv, nv, nt})
            for a = 0, nt do
                for b = 0, nv do
                    for c = b, nv do
                        q(b, c, a) = a + b - c + 1
                        q(c, b, a) = q(b, c, a)
                    end
                end
            end
            var f = [darray.DynamicMatrix(T)].new(&A, {nx, nv})
            f:fill([T](-1))
            var g = f:like()
            g:fill([T](2))
            
            var r = [darray.DynamicMatrix(T)].new(&A, {nx, nt})
            contraction.symcontract(&q, &f, &g, &r)

            var s = r:like()
            for l = 0, nx do
                for a = 0, nt do
                    s(l, a) = 0
                    for b = 0, nv do
                        for c = 0, nv do
                            s(l, a) = s(l, a) + q(b, c, a) * f(l, b) * g(l, c)
                        end -- for
                    end -- for
                end -- for
            end -- for
        end -- terracode
        
        test tmath.isapprox(&r, &s, tol)
    end
end
