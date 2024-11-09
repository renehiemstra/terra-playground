-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

import "terraform"
local err = require("assert")
local base = require("base")
local sarray = require("sarray")
local vecbase = require("vector")
local veccont = require("vector_contiguous")
local vecblas = require("vector_blas")
local concept = require("concept")
local sarray = require("sarray")

local StaticVector = terralib.memoize(function(T, N)
    
    --generate the raw type
    local SVector = sarray.SArrayRawType(T, {N})
    
    function SVector.metamethods.__typename(self)
        return ("StaticVector(%s, %d)"):format(tostring(T), N)
    end
    
    --add base functionality
    base.AbstractBase(SVector)

    --implement interfaces
    sarray.SArrayStackBase(SVector)
    sarray.SArrayVectorBase(SVector)
    sarray.SArrayIteratorBase(SVector)

    veccont.VectorContiguous:addimplementations{SVector}

    terra SVector:getblasinfo()
        return self:length(), &self.data, 1
    end

    vecblas.VectorBLAS:addimplementations{SVector}

    SVector.staticmethods.from = macro(
        function(...)
            local args = {...}
            assert(#args == N, "Length of input list does not match static dimension")
            local vec = symbol(SVector)
            local set_values = terralib.newlist()
            for i, v in ipairs(args) do
                set_values:insert(quote [vec].data[i - 1] = v end)
            end
            return quote
                var [vec] = SVector.new()
                [set_values]
            in
                [vec]
            end
        end
    )

    --vecbase.VectorBase(SVector)

    return SVector
end)

return {
    StaticVector = StaticVector
}
