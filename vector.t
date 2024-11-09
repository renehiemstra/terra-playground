-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

import "terraform"
local stack = require("stack")
local err = require("assert")
local template = require("template")
local concept = require("concept")
local tmath = require("mathfuns")

concept.Vector = concept.AbstractInterface:new("Vector")
local Stack = stack.Stack
concept.Vector:inheritfrom(Stack)
concept.Vector:addmethod{
	fill = concept.Number -> {},
	clear = {} -> {},
	sum  = {} -> concept.Number,
	-- BLAS operations
	copy = &Stack -> {},
	swap = &Stack -> {},
	scal = concept.Number -> {},
	axpy = {concept.Number, &Stack} -> {},
	dot = &Stack -> concept.Number,
	norm = {} -> concept.Number,
}

local VectorBase = function(Vector)

	assert(Stack(Vector),
		"A vector base implementation requires a valid stack implementation.")

    local T = Vector.eltype

    terra Vector:getbuffer()
        return self:length(), &self.data[0]
    end

    terra Vector:fill(value : T)
        for i = 0, self:length() do
            self:set(i, value)
        end
    end

    terraform Vector:copy(other : &S) where {S : Stack}
        for i = 0, self:length() do
            self:set(i, other:get(i))
        end
    end

	terraform Vector:swap(other : &S) where {S : Stack}
		err.assert(self:length() == other:length())
		for i = 0, self:length() do
			var tmp = other:get(i)
			other:set(i, self:get(i))
			self:set(i, tmp)
		end
	end

    if concept.Number(T) then

        terra Vector:clear()
		    self:fill(0)
	    end

        terra Vector:sum()
            var res : T = 0
            for i = 0, self:length() do
                res = res + self:get(i)
            end
            return res
        end

        terra Vector:scal(a : T)
            for i = 0, self:length() do
                self:set(i, a * self:get(i))
            end
        end

        terra Vector:axpy(a : T, x : &Vector)
            for i = 0, self:length() do
                self:set(i, self:get(i) + a * x:get(i))
            end
        end

        terra Vector:dot(other : &Vector)
            var s = T(0)
            for i = 0, self:length() do
                s = s + self:get(i) * other:get(i)
            end
            return s
        end
        
        terra Vector:norm2()
            return self:dot(self)
        end

        if concept.Float(T) then
            terra Vector:norm() : T
                return tmath.sqrt(self:norm())
            end
		end

    end

	if concept.Float(T) then
		assert(concept.Vector(Vector), "Incomplete implementation of vector base class")
		concept.Vector:addimplementations{Vector}
	end

end

return {
    Vector = concept.Vector,
    VectorBase = VectorBase
}
