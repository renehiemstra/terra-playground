-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

--lua function that generates a terra type that are function objects. these wrap
--a function in the 'apply' metamethod and store any captured variables in the struct
--as entries
local lambda_generator = function(fun, ...)
    --get the captured variables
    local captures = {...}
    --wrapper struct
    local lambda = terralib.types.newstruct("lambda")
    --add captured variable types as entries to the wrapper struct
    for i, sym in ipairs(captures) do
        lambda.entries:insert({field = "_"..tostring(i-1), type = sym.tree.type})
    end
    lambda:complete()
    --overloading the call operator - making 'lambda' a function object
    lambda.metamethods.__apply = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        local capt = terralib.newlist()
        for i,v in ipairs(self.tree.type.entries) do
            local field = "_"..tostring(i-1)
            capt:insert(quote in self.[field] end)
        end
        return `fun([args], [capt])
    end))
    --determine return-type from lambda expression
    lambda.returntype = fun.tree.type.type.returntype
    --determine parameter types and captured types
    local params = fun.tree.type.type.parameters
    local N, M = #params, #captures
    local K = N-M
    lambda.parameters, lambda.captures = terralib.newlist{}, terralib.newlist{}
    for k=1,K do
        lambda.parameters:insert(params[k])
    end
    for k=K+1,N do
        lambda.captures:insert(params[k])
    end
    --return function object
    return lambda
end

--return a function object with captured variables in ...
local lambda = macro(function(fun, ...)
    --get the captured variables
    local captures = {...}
    local p = lambda_generator(fun, ...)
    --create and return lambda object by value
    return quote
        var f = p{[captures]}
    in
        f
    end
end)

return {
    lambda = lambda,
    lambda_generator = lambda_generator 
}