local blas = require("axpy")
local mem = require("mem")
local complex = require("complex")
local C = terralib.includecstring[[
    #include <stdio.h>
]]
local complexFloat, If = unpack(complex(float))
local complexDouble, I = unpack(complex(double))

import "terratest/terratest"

for _, T in ipairs({float, double, complexFloat, complexDouble}) do
testenv(T) "BLAS level 1" do
        testset "swap scalar" do
            terracode
                var x = T(1)
                var y = T(2)
                blas.swap(1, &x, 1, &y, 1)
            end
            test x == T(2)
            test y == T(1)
        end -- testset

        testset "swap vectors" do
            local n = 4
            local incx = 3
            local incy = 2
            terracode
                var x: T[n * incx]
                var y: T[n * incy]
                var xold: T[n * incx]
                var yold: T[n * incy]
                for i = 0, n do
                    x[i * incx] = T(i)
                    xold[i * incx] = x[i * incx]
                    y[i * incy] = T(n - i)
                    yold[i * incy] = y[i * incy]
                end
                blas.swap(n, &x[0], incx, &y[0], incy)
            end

            for i = 0, n - 1 do
                test x[i * incx] == yold[i * incy]
                test y[i * incy] == xold[i * incx]
            end
        end -- testset swap

        testset "scal scalar" do
            terracode
                var x = T(2)
                var xold = x
                var a = T(6)
                blas.scal(1, a, &x, 1)
            end
            test x == xold * a
        end --testset scal

        testset "scal vector" do
            local n = 4
            local incx = 3
            terracode
                var x: T[n * incx]
                var xold: T[n * incx]
                var a = T(5)
                for i = 0, n do
                    x[i * incx] = T(i)
                    xold[i * incx] = x[i * incx]
                end
                blas.scal(n, a, &x[0], incx)
            end

            for i = 0, n - 1 do
                test x[i * incx] == a * xold[i * incx]
            end
        end -- testset scal

        testset "copy scalar" do
            terracode
                var x = T(3)
                var y = T(1)
                blas.copy(1, &x, 1, &y, 1)
            end

            test x == y
        end -- testset copy

        testset "copy vector" do
            local n = 4
            local incx = 3
            local incy = 2
            terracode
                var x: T[n * incx]
                var y: T[n * incy]

                for i = 0, n do
                    x[i * incx] = T(i)
                end
                blas.copy(n, &x[0], incx, &y[0], incy)
            end

            for i = 0, n - 1 do
                test y[i * incy] == x[i * incx]
            end
        end -- testset copy

        testset "axpy scalar" do
            terracode
                var a = T(2)
                var x = T(2)
                var y = T(3)
                var yold = y
                blas.axpy(1, a, &x, 1, &y, 1)
            end

            test y == a * x + yold
        end -- testset axpy

        testset "axpy vectors" do
            local n = 4
            local incx = 3
            local incy = 2
            terracode
                var a = T(2)
                var x: T[n * incx]
                var y: T[n * incy]
                var yold: T[n * incy]

                for i = 0, n do
                    x[i * incx] = T(i)
                    y[i * incy] = T(n - i)
                    yold[i * incy] = y[i * incy]
                end

                blas.axpy(n, a, &x[0], incx, &y[0], incy)
            end

            for i = 0, n - 1 do
                test y[i * incy] == a * x[i * incx] + yold[i * incy]
            end
        end -- testenv axpy

        testset "dot vectors real" do
            local n = 4
            local incx = 3
            local incy = 2
            terracode
                var x: T[n * incx]
                var y: T[n * incy]

                for i = 0, n do
                    x[i * incx] = T(i)
                    y[i * incy] = T(n - i)
                end
                var num = blas.dot(n, &x[0], incx, &y[0], incy)
                var ref = T(0)
                for i = 0, n do
                    ref = ref + x[i * incx] * y[i * incy]
                end
            end

            test ref == num
        end --testset dot
    end -- testenv BLAS level 1
end -- for
