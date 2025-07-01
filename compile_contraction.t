-- SPDX-FileCopyrightText: 2025 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2025 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

local compile = require("compile")
local contraction = require("contraction")
local darray = require("darray")

local terra bilinearcollisions(
    nx: int64,
    ntest: int64,
    ntrial: int64,
    qp: &double,
    fp: &double,
    gp: &double,
    rp: &double    
)
    var q = (
        [darray.DynamicArray(double, 3)]
            .frombuffer({ntrial, ntrial, ntest}, qp)
    )
    var f = [darray.DynamicMatrix(double)].frombuffer({nx, ntrial}, fp)
    var g = [darray.DynamicMatrix(double)].frombuffer({nx, ntrial}, gp)
    var r = [darray.DynamicMatrix(double)].frombuffer({nx, ntest}, rp)

    contraction.symcontract(&q, &f, &g, &r)
end

compile.generateCAPI(
    "contraction",
    {bilinearcollisions = bilinearcollisions}
)
