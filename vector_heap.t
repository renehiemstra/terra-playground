local alloc = require("alloc")
local err = require("assert")
local complex = require("complex")
local blas = require("blas")
local base = require("vector_base")
local stack = require("stack")

local complexFloat = complex.complex(float)
local complexDouble = complex.complex(double)

local function is_blas_type(T)
    local blas_type = {float, double, complexFloat, complexDouble}
    for _, B in pairs(blas_type) do
        if B == T then
            return true
        end
    end

    return false
end

local VectorHeap = terralib.memoize(function(T, A)
    A = A or alloc.Default
    alloc.Allocater:isimplemented(A)
	local S = stack.DynamicStack(T, int64, A)

    local struct vector{
		data: S
        inc: int64
    }

    local terra new(size: int64)
		return vector {S.new(size), 1}
    end

    terra vector:free()
		self.data:free()
    end

    terra vector:size()
        return self.data:size()
    end

    terra vector:get(i: int64)
        err.assert(i >= 0)
        err.assert(i < self:size())
		return self.data:get(self.inc * i)
    end

    terra vector:set(i: int64, a: T)
        err.assert(i >= 0)
        err.assert(i < self:size())
		self.data:set(self.inc * i, a)
    end

    vector = base.VectorBase(vector, T, int64)

	terra vector:push(a: T)
		self.data:push(a)
	end

	terra vector:pop()
		return self.data:pop()
	end

    terra vector:subview(size: int64, offset: int64, inc: int64)
		return vector {self.data:slice(size, offset * self.inc), inc * self.inc}
    end

    terra vector:inc()
        return self.inc
    end

    terra vector:data()
        return self.data:data()
    end

    terra vector:getblasinfo()
        return self:size(), self:data(), self:inc()
    end

	--[[
    if is_blas_type(T) then
        terra vector:swap(x: &vector)
            var x_size, x_data, x_inc = x:getblasinfo()
            var y_size, y_data, y_inc = self:getblasinfo()

            err.assert(x_size == y_size)

            blas.swap(x_size, x_data, x_inc, y_data, y_inc)
        end

        terra vector:scal(a: T)
            var x_size, x_data, x_inc = self:getblasinfo()

            blas.scal(x_size, a, x_data, x_inc)
        end
        
        terra vector:axpy(a: T, x: &vector)
            var x_size, x_data, x_inc = x:getblasinfo()
            var y_size, y_data, y_inc = self:getblasinfo()

            err.assert(x_size == y_size)

            blas.axpy(x_size, a, x_data, x_inc, y_data, y_inc)
        end

        terra vector:dot(x: &vector)
            var x_size, x_data, x_inc = x:getblasinfo()
            var y_size, y_data, y_inc = self:getblasinfo()

            err.assert(x_size == y_size)

            return blas.dot(x_size, x_data, x_inc, y_data, y_inc)
        end

        terra vector:nrm2()
            var x_size, x_data, x_inc = self:getblasinfo()

            return blas.nrm2(x_size, x_data, x_inc)
        end

        terra vector:asum()
            var x_size, x_data, x_inc = self:getblasinfo()

            return blas.asum(x_size, x_data, x_inc)
        end

        terra vector:iamax()
            var x_size, x_data, x_inc = self:getblasinfo()

            return blas.iamax(x_size, x_data, x_inc)
        end
    end
	--]]

    vector.metamethods.__for = function(iter, body)
        return quote
            var size = iter:size()
            for i = 0, size do
                var data = iter:get(i)
                [body(data)]
            end
        end
    end

    local from = macro(
        function(...)
            local arg = {...}
            local vec = symbol(vector)
            local push = terralib.newlist()
            for _, v in ipairs(arg) do
                push:insert(quote [vec]:push(v) end)
            end

            return quote
                       var [vec] = new(0)
                       [push]     
                   in
                       [vec]
                   end
        end)

    local terra from_buffer(size: int64, data: &T, inc: int64)
		return vector {S.frombuffer(size * inc, data), inc}
    end

    local terra like(x: vector)
        var size = x:size()
        var y = new(size)

        return y
    end

    local terra zeros_like(x: vector)
        var y = like(x)
        y:clear()

        return y
    end

    local static_methods = {
        new = new,
        from = from,
        frombuffer = from_buffer,
        like = like,
        zeroslike = zeros_like
	}

	-- TODO is this the right place? Or should this be in a meta type?
	vector.metamethods.__methodmissing = macro(function(name, obj, ...)
		local is_static = (static_methods[name] ~= nil)
		local args = terralib.newlist({...})
		if is_static then
			args:insert(1, obj)
		end
		local types = args:map(function(t) return t.tree.type end)
		if is_static then
			local method = static_methods[name]
			return `method([args])
		else
			types:insert(1, &vector)
			local method = vector.template[name]
			print("Types are", types)
			local func = method(unpack(types))
			return quote var self = obj in [func](&self, [args]) end
		end
	end)

    return vector
end)

return {
    VectorHeap = VectorHeap
}
