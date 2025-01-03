-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

import "terratest/terratest"

local lu = require("lu")
local alloc = require("alloc")
local concepts = require("concepts")
local random = require("random")
local complex = require("complex")
local nfloat = require("nfloat")
local dvector = require("dvector")
local dmatrix = require("dmatrix")
local mathfun = require("mathfuns")

local float128 = nfloat.FixedFloat(128)
local float1024 = nfloat.FixedFloat(1024)
local tol = {["float"] = 1e-6,
             ["double"] = 1e-15,
             [tostring(float128)] = `[float128]("1e-30"),
             [tostring(float1024)] = `[float1024]("1e-300"),
            }

for _, Ts in pairs({float, double, float128, float1024}) do
    for _, is_complex in pairs({false, true}) do
        local T = is_complex and complex.complex(Ts) or Ts
        local unit = is_complex and T:unit() or 0
        local DMat = dmatrix.DynamicMatrix(T)
        local DVec = dvector.DynamicVector(T)
        local PVec = dvector.DynamicVector(int32)
        local Alloc = alloc.DefaultAllocator(Ts)
        local Rand = random.Default(float)
        local LUDense = lu.LUFactory(DMat, PVec)

        testenv(T) "LU factorization for small matrix" do
            terracode
                var alloc: Alloc
                var a = DMat.from(&alloc, {{1, 2}, {3, 4}})
                var p = PVec.zeros(&alloc, a:rows())
                var tol: Ts = [ tol[tostring(Ts)] ]
                var lu = LUDense.new(&a, &p, tol)
                var x = DVec.from(&alloc, 2, 1)
                var xt = DVec.from(&alloc, 4, 2)
                lu:factorize()
                lu:solve(false, &x)
                lu:solve(true, &xt)
            end

            if not concepts.BLASNumber(T) then
                testset "Factorize" do
                    test mathfun.abs(a(0, 0) - [T](3)) < tol
                    test mathfun.abs(a(0, 1) - [T](4)) < tol
                    test mathfun.abs(a(1, 0) - [T](1) / 3) < tol
                    test mathfun.abs(a(1, 1) - [T](2) / 3) < tol

                    test p(0) == 1
                    test p(1) == 0
                end
            end

            testset "Solve" do
                test mathfun.abs(x(0) - [T](-3)) < tol
                test mathfun.abs(x(1) - [T](5) / 2) < tol
            end

            testset "Solve transposed" do
                test mathfun.abs(xt(0) - [T](-5)) < tol
                test mathfun.abs(xt(1) - [T](3)) < tol
            end
        end

        testenv(T) "LU factorization for random matrix" do
            local n = 41
            terracode
                var alloc: Alloc
                var rand = Rand.from(2359586)
                var a = DMat.new(&alloc, n, n)
                var x = DVec.new(&alloc, n)
                var y = DVec.zeros_like(&alloc, &x)
                var yt = DVec.zeros_like(&alloc, &x)
                for i = 0, n do
                    x(i) = rand:rand_normal(0, 1) + [unit] * rand:rand_normal(0, 1)
                    for j = 0, n do
                        a(i, j) = rand:rand_normal(0, 1) + [unit] * rand:rand_normal(0, 1)
                    end
                end
                a:apply(false, [T](1), &x, [T](0), &y)
                a:apply(true, [T](1), &x, [T](0), &yt)
                var p = PVec.zeros(&alloc, n)
                var tol: Ts = [ tol[tostring(Ts)] ]
                var lu = LUDense.new(&a, &p, tol)
                lu:factorize()
                lu:solve(false, &y)
                lu:solve(true, &yt)
            end

            testset "Solve" do
                for i = 0, n - 1 do
                    test mathfun.abs(y(i) - x(i)) < 1000 * tol * mathfun.abs(x(i)) + tol
                end
            end

            testset "Solve transposed" do
                for i = 0, n - 1 do
                    test mathfun.abs(yt(i) - x(i)) < 2000 * tol * mathfun.abs(x(i)) + tol
                end
            end
        end
    end
end
