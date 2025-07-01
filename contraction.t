-- SPDX-FileCopyrightText: 2025 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2025 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT
local alloc = require("alloc")
local concepts = require("concepts")
local lambda = require("lambda")
local range = require("range")
local simd = require("simd")
local thread = require("thread")

import "terraform"

local Unit = range.Unitrange(uint64)
local symcontract
terraform symcontract(q: &Q, f: &F, g: &G, r: &R) where {
    Q: concepts.Tensor(concepts.Number, 3),
    F: concepts.Matrix(concepts.Number),
    G: concepts.Matrix(concepts.Number),
    R: concepts.Matrix(concepts.Number)
}
    var rn = Unit.new(0, f:rows())
    var A: alloc.DefaultAllocator()
    var go = lambda.new(
        [
            terra(l: uint64, q: q.type, f: f.type, g: g.type, r: r.type)
                var nt = q:size()[2]
                var nv = q:size()[0]
                for a = 0, nt do
                    var rr = [R.traits.eltype](0)
                    for b = 0, nv do
                        rr = rr + q(b, b, a) * f(l, b) * g(l, b)
                        for c = b + 1, nv do
                            rr = rr + 2 * q(b, c, a) * f(l, b) * g(l, c)
                        end
                    end
                    r(l, a) = rr
                end
            end
        ],
        {q = q, f = f, g = g, r = r}
    )
    thread.parfor(&A, rn, go)
end

do
    local function single_contract(T, N)
        local SIMD = vector(T, N)
        local terra lowlevelcontract(
            ntest: int64,
            ntrial: int64,
            q: &T,
            offset: int64,
            f: &T,
            g: &T,
            r: &T
        )
            var rr: SIMD = [T](0)
            for b = 0, ntrial do
                var idx = ntrial * ntest * b + ntest * b
                var qq = [simd.load(T, N)](q + offset + idx)
                rr = rr + qq * f[b] * g[b]
                for c = b + 1, ntrial do
                    var idx = ntrial * ntest * b + ntest * c
                    var qq = [simd.load(T, N)](q + offset + idx)
                    rr = rr + [T](2) * qq * f[b] * g[c]
                end
            end
            [simd.store(T, N)](r + offset, rr)
        end
        return lowlevelcontract
    end

    local C = terralib.includec("stdio.h")
    local MAX_POWER = 7
    terraform symcontract(q: &Q, f: &F, g: &G, r: &R) where {
        Q: concepts.Tensor(concepts.BLASFloat, 3),
        F: concepts.Matrix(concepts.BLASFloat),
        G: concepts.Matrix(concepts.BLASFloat),
        R: concepts.Matrix(concepts.BLASFloat)
    }
        var go = lambda.new([
            terra(l: uint64, q: q.type, f: f.type, g: g.type, r: r.type)
                var ntest = q:size()[2]
                var ntrial = q:size()[0]
                var start: int64 = 0
                var vecsize: int64 = [2 ^ MAX_POWER]
                var size: int64 = ntest
                while size ~= 0 do
                    while vecsize > size do
                        vecsize = vecsize / 2
                    end
                    escape
                        local T = Q.traits.eltype
                        for i = 0, MAX_POWER do
                            local N = 2 ^ i
                            local func = single_contract(T, N)
                            emit quote
                                if vecsize == N then
                                    func(
                                        ntest,
                                        ntrial,
                                        &q(0, 0, 0),
                                        start,
                                        &f(l, 0),
                                        &g(l, 0),
                                        &r(l, 0)
                                    )
                                end -- if
                            end -- quote
                        end -- for
                    end -- escape
                    size = size - vecsize
                    start = start + vecsize
                end -- while
            end -- terra
            ],
            {q = q, f = f, g = g, r = r}
        )
        var A: alloc.DefaultAllocator()
        var rn = Unit.new(0, f:rows())
        thread.parfor(&A, rn, go)
    end
end

return {
    symcontract = symcontract
}

