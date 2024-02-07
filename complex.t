local sqrt = terralib.overloadedfunction("sqrt")
for tname, ttype in pairs{f32 = float, f64 = double} do
    local d = terra(x: ttype) return [terralib.intrinsic("llvm.sqrt."..tname, ttype -> ttype)](x) end
    sqrt:adddefinition(d)
end

local function complex(T)

    local struct complex{
        re: T
        im: T
    }
    complex:setconvertible("array")

    terra complex.metamethods.__add(self: complex, other: complex)
        return complex {self.re + other.re, self.im + other.im}
    end

    terra complex.metamethods.__mul(self: complex, other: complex)
        return complex {self.re * other.re - self.im * other.im, self.re * other.im + self.im * other.re}
    end

    terra complex.metamethods.__unm(self: complex)
        return complex {-self.re, -self.im}
    end

    terra complex:normsq(): T
        return self.re * self.re + self.im * self.im
    end

    terra complex:norm(): T
        return sqrt(self:normsq())
    end

    terra complex:inverse()
       var nrmsq = self:normsq()
       return {self.re / nrmsq, -self.im / nrmsq}
    end

    terra complex.metamethods.__sub(self: complex, other: complex)
        return self + (-other)
    end

    terra complex.metamethods.__div(self: complex, other: complex)
        return self * other:inverse()
    end


    return complex
end

complex = terralib.memoize(complex)


return complex
