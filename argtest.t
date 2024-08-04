local io = terralib.includec("stdio.h")


local struct interval{
    a : int
    b : int
}

local min = terra(a : int, b : int)
    return terralib.select(a < b, a, b)
end

local max = terra(a : int, b : int)
    return terralib.select(a > b, a, b)
end

local intersect = terra(I : interval, J : interval)
    var a = max(I.a, J.a)
    var b = min(I.b, J.b)
    return interval{a, b}    
end


terra main()

    var I = interval{0,1}
    var J = interval{1,2}
    var K = intersect(I,J)
    io.printf("{a, b} = {%d, %d}\n", K.a, K.b)
end
main()