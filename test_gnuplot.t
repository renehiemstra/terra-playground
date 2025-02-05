-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

import "terratest/terratest"

local gnuplot = require("gnuplot")


if not __silent__ then

    terra test1()
        var fig : gnuplot.handle
        --gnuplot.setterm(fig, "wxt", 600, 400)
        gnuplot.plot_equation(fig, "sin(x)", "Sine wave")
    end
    test1()

    terra test2()
        var fig : gnuplot.handle
        --gnuplot.setterm(fig, "wxt", 600, 400)
        gnuplot.cmd(fig, "set samp 1000")
        gnuplot.cmd(fig, "set xrange [-.05:20]")
        gnuplot.cmd(fig, "set y2tics nomirror")
        gnuplot.cmd(fig, "set ytics nomirror")
        gnuplot.cmd(fig, "set ylabel 'Y0'")        
        gnuplot.cmd(fig, "set y2label 'J0'")
        gnuplot.cmd(fig, "set grid")
        gnuplot.cmd(fig, "plot besy0(x) axes x1y1 lw 2 title 'Y0', besj0(x) axes x1y2 lw 2 title 'J0'")
    end
    test2()

--[[
    terra test3()
        var fig : gnuplot.handle

        var x : double[50]
        var y : double[50]
        for i = 0, 50 do
            x[i] = i / 10.0
            y[i] = x[i] * x[i]
        end
        gnuplot.plot_coordinates(fig.h, &x[0], &y[0], 50, "parabola")

    end
    test3()
--]]
end