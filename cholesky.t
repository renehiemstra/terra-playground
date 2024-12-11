-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

import "terraform"

local base = require("base")
local err = require("assert")
local concepts = require("concepts")
local matbase = require("matrix")
local vecbase = require("vector")
local matblas = require("matrix_blas_dense")
local mathfun = require("mathfuns")
local lapack = require("lapack")

local Bool   = concepts.Bool
local Number = concepts.Number
local Vector = concepts.Vector(Number)
local Matrix = concepts.Matrix(Number)

local BLASNumber = concepts.BLASNumber
local BLASVector = concepts.BLASVector(BLASNumber)
local BLASMatrix = concepts.BLASDenseMatrix(BLASNumber)

terraform factorize(a: &M, tol: T) where {M: Matrix, T: Number}
    var n = a:rows()
    for i = 0, n do
        for j = 0, i + 1 do
            var sum = a:get(i, j)
            for k = 0, j do
                sum = sum - a:get(i, k) * mathfun.conj(a:get(j, k))
            end
            if i == j then
                var sumabs = mathfun.abs(sum)
                err.assert(mathfun.abs(sum - sumabs) < tol * sumabs + tol)
                a:set(i, i, mathfun.sqrt(sumabs))
            else
                a:set(i, j, sum / a:get(j, j))
            end
        end
    end
end

terraform factorize(a: &M, tol: T) where {M: BLASMatrix, T: BLASNumber}
    var n, m, adata, lda = a:getblasdenseinfo()
    err.assert(n == m)
    lapack.potrf(lapack.ROW_MAJOR, @"L", n, adata, lda)
end

local conj = mathfun.conj
terraform solve(trans: B, a: &M, x: &V) where {B: Bool, M: Matrix, V: Vector}
    var n = a:rows()
    for i = 0, n do
        for k = 0, i do
            x:set(i, x:get(i) - a:get(i, k) * x:get(k))
        end
        x:set(i, x:get(i) / a:get(i, i))
    end

    for ii = 0, n do
        var i = n - 1 - ii
        for k = i + 1, n do
            x:set(i, x:get(i) - mathfun.conj(a:get(k, i)) * x:get(k))
        end
        x:set(i, x:get(i) / a:get(i, i))
    end
end

terraform solve(trans: B, a: &M, x: &V)
    where {B: Bool, M: BLASMatrix, V: BLASVector}
    var n, m, adata, lda = a:getblasdenseinfo()
    err.assert(n == m)
    var nx, xdata, incx = x:getblasinfo()
    lapack.potrs(lapack.ROW_MAJOR, @"L", n, 1, adata, lda, xdata, incx)
end

local CholeskyFactory = terralib.memoize(function(M)

    local T = M.eltype
    local Vector = concepts.Vector(T)
    local Matrix = concepts.Matrix(T)
    local Factorization = concepts.Factorization(T)

    assert(Matrix(M), "Type " .. tostring(M)
                              .. " does not implement the matrix interface")
    local Ts = T
    local Ts = concepts.Complex(T) and T.traits.eltype or T
    local struct cho{
        a: &M
        tol: Ts
    }
    function cho.metamethods.__typename(self)
        return ("CholeskyFactorization(%s)"):format(tostring(T))
    end
    base.AbstractBase(cho)

    terra cho:rows()
        return self.a:rows()
    end

    terra cho:cols()
        return self.a:cols()
    end

    terra cho:factorize()
        escape
            local impl = factorize:dispatch(&M, Ts)
            emit quote return impl(self.a, self.tol) end
        end
    end

    terraform cho:solve(trans: bool, x: &V) where {V: Vector}
        return solve(trans, self.a, x)
    end

    terraform cho:apply(trans: B, a: T, x: &V1, b: T, y: &V2)
        where {B: Bool, V1: Vector, V2: Vector}
        self:solve(trans, x)
        y:scal(b)
        y:axpy(a, x)
    end

    assert(Factorization(cho))

    cho.staticmethods.new = terra(a: &M, tol: Ts)
        err.assert(a:rows() == a:cols())
        return cho {a, tol}
    end

    return cho
end)

return {
    CholeskyFactory = CholeskyFactory,
}
