-- Wrap FLINT without inlines
local flint = terralib.includec("flint/nfloat.h", {"-DNFLOAT_INLINES_C=1"})
local gr = terralib.includec("flint/gr.h", {"-DGR_INLINES_C=1"})
terralib.linklibrary("libflint.so")
local io = terralib.includec("stdio.h")
local mathfun = require("mathfuns")

local suffix = {64, 128, 192, 256, 384, 512, 1024, 2048, 4096}
local float_type = {}
local context = {}
for _, N in pairs(suffix) do
    float_type[N] = flint[string.format("nfloat%d_struct", N)]
    -- Meta information on fixed precision floats is stored in a context.
    -- Mathematically, they represent rings.
    -- Here, we store them as global variables in table such that
    -- each float type has exactly one context it will use.
    context[N] = global(flint.gr_ctx_t)
    local ctx = context[N]:get()
    -- Call clean_context() to release memory allocated by nfloat_ctx_init
    flint.nfloat_ctx_init(ctx, N, 0)
end

local unary_math = {
    "abs",
    "sqrt",
    "floor",
    "ceil",
    "exp",
    "expm1",
    "log",
    "log1p",
    "sin",
    "cos",
    "tan",
    "sinh",
    "cosh",
    "tanh",
    "gamma",
}

local binary_math = {
    "pow",
}

local FixedFloat = terralib.memoize(function(N)
    local ctype = float_type[N]
    assert(ctype, "No support for precision " .. N .. " in FixedFloat")
    local ctx = context[N]:get()

    local struct nfloat {
        data: ctype
    }

    function nfloat.metamethods.__typename()
        return string.format("FixedFloat(%d)", N)
    end

    local terra new()
        var data: ctype
        flint.nfloat_init(&data, ctx)
        return nfloat {data}
    end

    local terra from_double(x: double)
        var f = new()
        flint.nfloat_set_d(&f.data, x, ctx)
        return f
    end
    local terra from_str(s: rawstring)
        var f = new()
        flint.nfloat_set_str(&f.data, s, ctx)
        return f
    end
    local from = terralib.overloadedfunction("from", {from_double, from_str})

    local to_str = macro(function(x)
        local digits = math.floor(N * (math.log(2) / math.log(10)))
        return quote
                var str: rawstring
                -- TODO: Fix memory leak
                -- defer flint.flint_free(str)
                gr.gr_get_str_n(&str, &x.data, [digits], ctx)
            in
                str
        end        
    end)

    function nfloat.metamethods.__cast(from, to, exp)
        if to == nfloat then
            if from:isarithmetic() then
                return `from_double(exp)
            elseif from:ispointer() and from.type == int8 then
                return `from_str(exp)
            else
                error("Cannot cast from " .. from .. " to " .. to)
            end
        end
        error("Unknown type")
    end

    local binary = {
        __add = flint.nfloat_add,
        __mul = flint.nfloat_mul,
        __sub = flint.nfloat_sub,
        __div = flint.nfloat_div,
    }
    for key, method in pairs(binary) do
        nfloat.metamethods[key] = terra(self: nfloat, other:nfloat)
            var res = new()
            [method](&res.data, &self.data, &other.data, ctx)
            return res
        end
    end

    local unary = {
        __unm = flint.nfloat_neg,
    }
    for key, method in pairs(unary) do
        nfloat.metamethods[key] = terra(self: nfloat)
            var res = new()
            [method](&res.data, &self.data, ctx)
            return res
        end
    end

    local boolean = {
        __eq = flint.nfloat_equal
    }

    for key, method in pairs(boolean) do
        nfloat.metamethods[key] = terra(self: nfloat, other: nfloat)
            return [method](&self.data, &other.data, ctx)
        end
    end

    local terra pi()
        var res = new()
        flint.nfloat_pi(&res, ctx)
        return res
    end

    for _, func in pairs(unary_math) do
        local name = "nfloat_" .. func
        local terra impl(x: nfloat)
            var y: nfloat
            flint.[name](&y.data, &x.data, ctx)
            return y
        end
        mathfun[func]:adddefinition(impl)
    end

    for _, func in pairs(binary_math) do
        local name = "nfloat_" .. func
        local terra impl(x: nfloat, y: nfloat)
            var z: nfloat
            flint.[name](&z.data, &x.data, &y.data, ctx)
            return z
        end
        mathfun[func]:adddefinition(impl)
    end

    mathfun.fusedmuladd:adddefinition(
        terra impl(x: nfloat, y: nfloat, z: nfloat)
            return x * y + z
        end
    )

    local staticmethods = {
        from = from,
        tostr = to_str,
        pi = pi,
    }

    nfloat.metamethods.__getmethod = function(self, methodname)
        return staticmethods[methodname] or nfloat.methods[methodname]
    end

    return nfloat
end)

local terra clean_context()
    escape
        for _, N in pairs(suffix) do
            local val = context[N]:get()
            emit quote
                gr.gr_ctx_clear(val)
            end
        end
    end
end

return {
    FixedFloat = FixedFloat,
    clean_context = clean_context
}
