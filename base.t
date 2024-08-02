local concept = require("concept")

local Base = {}

function Base:new(name, func)
	local base = {name = name}
	setmetatable(base, {__index = self})

	local mt = getmetatable(base)

	function mt:__call(T)
		func(T)
	end

	function mt.__mul(B1, B2)
		local function impl(T)
			B1(T)
			B2(T)
		end
		return Base:new(B1.name .. "And" .. B2.name, impl)
	end

	return base
end

local AbstractBase = Base:new("AbstractBase",
	function(T)
		assert(terralib.types.istype(T))
		assert(T:isstruct())
		local Self = concept.Concept:new(tostring(T),
										 function(Tp) print(T.name, Tp.name, Tp.name == T.name); return Tp.name == T.name end
										)
		local SelfPtr = concept.Ptr(Self)
		for key, val in pairs({staticmethods = {}, templates = {}, Self = Self,
													 SelfPtr = SelfPtr})  do
			if T.key == nil then
				rawset(T, key, val)
			end
		end
		Self.methods = T.methods
		Self.staticmethods = T.staticmethods
		Self.templates = T.templates

		T.metamethods.__getmethod = function(self,methodname)
	    local fnlike = self.methods[methodname] or self.staticmethods[methodname]
	    if not fnlike and terralib.ismacro(self.metamethods.__methodmissing) then
	        fnlike = terralib.internalmacro(function(ctx,tree,...)
	            return self.metamethods.__methodmissing:run(ctx,tree,methodname,...)
	        end)
	    end
	    return fnlike
		end

		T.metamethods.__methodmissing = macro(function(methodname,...)
    local gen = X.template[methodname]
    if gen then
        local args = terralib.newlist{...}
        local types = args:map(function(v) return v.tree.type end)
        local f = gen(unpack(types))
        return `f(args)
    end
end)

		T.metamethods.__methodmissing = macro(function(name, obj, ...)
			local is_static = (T.staticmethods[name] ~= nil)
			local args = terralib.newlist({...})
			local types = args:map(function(t) return t.tree.type end)
			if is_static then
				local method = T.staticmethods[name]
				return quote method([args]) end
			else
				types:insert(1, &T)
				local method = T.templates[name]
				local func = method(unpack(types))
				return quote [func](&obj, [args]) end
			end
		end)
	end
)

return {
	Base = Base,
	AbstractBase = AbstractBase
}
