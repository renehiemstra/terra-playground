


local interval = Interval(double)
local hypercube = Hypercube(double, 2, 1)

terra main()
    var cube_1 = hypercube.new(interval.new(1,2), 0)
    var cube_2 = hypercube.new(0, interval.new(1,2))
end
main()


testenv "pyramid - N = 2" do

    local point = squad.Point(double, 2)
    local interval = squad.Interval(double)
    local hypercube = squad.Hypercube(double, 2)
    local pyramid = squad.Pyramid(double, 2)
    local pyramids = squad.Pyramids(double, 2)
    
    terracode
        var e1 = hypercube.new(interval.new(1, 3), interval.new(5, 5))
        var e2 = hypercube.new(interval.new(3, 3), interval.new(2, 5))
        var apex = point.from(1,2)
        var pyr_a = pyramid.new(e2, apex)
        var pyr_b = pyramid.new(e1, apex)
    end
    test e1:vol()==2 and e2:vol()==3

    testset "new, dim, vol, height" do
        --first pyramid
        test pyr_a:dim()==2
        test pyr_a:height()==2
        --test pyr_a:vol()==0.5 * e1:vol() * e2:vol()
        --second pyramid
        --test pyr_b:dim()==2
        --test pyr_b:height()==3
        --test pyr_b:vol()==0.5 * e1:vol() * e2:vol()
    end

    testset "evaluation" do
        --first pyramid
        test pyr_a(0,0) == apex
        test pyr_a(0,1) == point.from(3,2)
        test pyr_a(1,1) == point.from(3,5)
        --second pyramid
        test pyr_b(0,0) == apex
        test pyr_b(0,1) == point.from(1,5)
        test pyr_b(1,1) == point.from(3,5)
    end

    testset "jacobian" do
        --first pyramid
        test pyr_a:jac(0,0) == 0
        test pyr_a:jac(0,1) == 6
        test pyr_a:jac(1,1) == 6
        --second pyramid
        test pyr_b:jac(0,0) == 0
        test pyr_b:jac(0,1) == 6
        test pyr_b:jac(1,1) == 6
    end

    testset "decompose square into pyramids" do
        terracode
            var count = 0
            for pyr in pyramids{e1*e2, apex} do
                count = count + 1
            end
        end
        test count == 2
    end

end

testenv "pyramid - N = 3" do

    local point = squad.Point(double, 3)
    local interval = squad.Interval(double)
    local hypercube = squad.Hypercube(double, 3)
    local pyramid = squad.Pyramid(double, 3)
    local pyramids = squad.Pyramids(double, 3)
    
    terracode
        var base = hypercube.new(interval.new(1, 2), interval.new(1, 3), interval.new(1, 1))
        var apex = point.from(0,0,2)
        var pyr = pyramid.new(base, apex)
    end
    test base:vol()==2

    testset "new, dim, vol, height" do
        test pyr:dim()==3
        test pyr:height()==1
        test pyr:vol()==2. / 3.
    end

    testset "evaluation" do
        test pyr(0,0,1)==base(0,0)
        test pyr(0,0,0)==apex
        test pyr(0.0,0.0,0.5) == point.from(0.5,0.5,1.5)
        test pyr(0.5,0.5,0.5) == point.from(3./4.,1,1.5)
    end

    testset "jacobian" do
        test pyr:jac(0,0,1) == 2
        test pyr:jac(0.0,0.0,0.5) == 0.5
        test pyr:jac(0.5,0.5,0.5) == 0.5
        test pyr:jac(0,0,0) == 0
    end

    testset "decompose square into pyramids" do
        terracode
            var cube = hypercube.new(interval.new(1, 2), interval.new(1, 3), interval.new(1, 2))
            var apex = point.from(1,1,1)
            var count = 0
            var volume = 0.0
            for pyr in pyramids{cube, apex} do
                volume = volume + pyr:vol()
                count = count + 1
            end
        end
        test count == 3
        test volume == 2
    end

end

testenv "product-pair - N = 3, D = 0" do

    local point = squad.Point(double, 3)
    local interval = squad.Interval(double)
    local hypercube = squad.Hypercube(double, 3)
    local product = squad.ProductPair(double, 3)
    
    terracode
        var A = hypercube.new(interval.new(0, 1), interval.new(0, 1), interval.new(0, 1))
        var B = hypercube.new(interval.new(1, 2), interval.new(1, 2), interval.new(1, 2))
        var C = hypercube.intersect(A, B)
        var X = product.new(C, A / C, point.from(1,1,1))
        var Y = product.new(C, B / C, point.from(1,1,1))
    end

    testset "new, vol, origin" do
        --X
        test X.a:vol()==1 and X.b:vol()==1
        test X.a:getorigin() == point.from(1,1,1)
        test X.b:getorigin() == point.from(1,1,1)
        --Y
        test Y.a:vol()==1 and Y.b:vol()==1
        test Y.a:getorigin() == point.from(1,1,1)
        test Y.b:getorigin() == point.from(1,1,1)
    end

    testset "evaluation" do
        --X
        test X({},{0,0,0}) == X:getorigin()
        test X({},{1,1,1}) == point.from(0,0,0)
        --Y
        test Y({},{0,0,0}) == Y:getorigin()
        test Y({},{1,1,1}) == point.from(2,2,2)
    end
end

testenv "product-pair - N = 3, D = 1" do

    local point = squad.Point(double, 3)
    local interval = squad.Interval(double)
    local hypercube = squad.Hypercube(double, 3)
    local product = squad.ProductPair(double, 3)
    
    terracode
        var A = hypercube.new(interval.new(0, 1), interval.new(0, 1), interval.new(0, 1))
        var B = hypercube.new(interval.new(1, 2), interval.new(1, 2), interval.new(0, 1))
        var C = hypercube.intersect(A, B)
        var X = product.new(C, A / C, point.from(1,1,0))
        var Y = product.new(C, B / C, point.from(1,1,0))
    end

    testset "new, vol, origin" do
        --X
        test X.a:vol()==1 and X.b:vol()==1
        test X.a:getorigin() == point.from(1,1,0)
        test X.b:getorigin() == point.from(1,1,0)
        --Y
        test Y.a:vol()==1 and Y.b:vol()==1
        test Y.a:getorigin() == point.from(1,1,0)
        test Y.b:getorigin() == point.from(1,1,0)
    end

    testset "evaluation" do
        --X
        test X({0},{0,0}) == X:getorigin()
        test X({0},{1,0}) == point.from(0, 1, 0)
        test X({0},{0,1}) == point.from(1, 0, 0)
        test X({0},{1,1}) == point.from(0, 0, 0)
        test X({1},{0,0}) == point.from(1, 1, 1)
        test X({1},{1,0}) == point.from(0, 1, 1)
        test X({1},{0,1}) == point.from(1, 0, 1)
        test X({1},{1,1}) == point.from(0, 0, 1)
        --Y
        test Y({0},{0,0}) == Y:getorigin()
        test Y({0},{1,0}) == point.from(2, 1, 0)
        test Y({0},{0,1}) == point.from(1, 2, 0)
        test Y({0},{1,1}) == point.from(2, 2, 0)
        test Y({1},{0,0}) == point.from(1, 1, 1)
        test Y({1},{1,0}) == point.from(2, 1, 1)
        test Y({1},{0,1}) == point.from(1, 2, 1)
        test Y({1},{1,1}) == point.from(2, 2, 1)
    end
end

testenv "product-pair - N = 3, D = 2" do

    local point = squad.Point(double, 3)
    local interval = squad.Interval(double)
    local hypercube = squad.Hypercube(double, 3)
    local product = squad.ProductPair(double, 3)
    
    terracode
        var A = hypercube.new(interval.new(1, 3), interval.new(0, 1), interval.new(0, 1))
        var B = hypercube.new(interval.new(3, 4), interval.new(0, 1), interval.new(0, 1))
        var C = hypercube.intersect(A, B)
        var X = product.new(C, A / C, point.from(3,0,0))
        var Y = product.new(C, B / C, point.from(3,0,0))
    end

    testset "new, vol, origin" do
        --X
        test X.a:vol()==1 and X.b:vol()==2
        test X.a:getorigin() == point.from(3,0,0)
        test X.b:getorigin() == point.from(3,0,0)
        --Y
        test Y.a:vol()==1 and Y.b:vol()==1
        test Y.a:getorigin() == point.from(3,0,0)
        test Y.b:getorigin() == point.from(3,0,0)
    end

    testset "evaluation" do
        --X
        test X({0,0},{0}) == X:getorigin()
        test X({1,0},{0}) == point.from(3, 1, 0)
        test X({0,1},{0}) == point.from(3, 0, 1)
        test X({1,1},{0}) == point.from(3, 1, 1)
        test X({1,0},{1}) == point.from(1, 1, 0)
        test X({0,1},{1}) == point.from(1, 0, 1)
        test X({1,1},{1}) == point.from(1, 1, 1)
        --Y
        test Y({0,0},{0}) == Y:getorigin()
        test Y({1,0},{0}) == point.from(3, 1, 0)
        test Y({0,1},{0}) == point.from(3, 0, 1)
        test Y({1,1},{0}) == point.from(3, 1, 1)
        test Y({1,0},{1}) == point.from(4, 1, 0)
        test Y({0,1},{1}) == point.from(4, 0, 1)
        test Y({1,1},{1}) == point.from(4, 1, 1)
    end
end

testenv "product-pair - N = 3, D = 3" do

    local point = squad.Point(double, 3)
    local interval = squad.Interval(double)
    local hypercube = squad.Hypercube(double, 3)
    local product = squad.ProductPair(double, 3)
    
    terracode
        var A = hypercube.new(interval.new(0, 1), interval.new(0, 1), interval.new(0, 1))
        var B = hypercube.new(interval.new(0, 1), interval.new(0, 1), interval.new(0, 1))
        var C = hypercube.intersect(A, B)
        var X = product.new(C, A / C, point.from(0,0,0))
        var Y = product.new(C, B / C, point.from(0,0,0))
    end

    testset "new, vol, origin" do
        --X
        test X.a:vol()==1 and X.b:vol()==1
        test Y.a:getorigin() == point.from(0,0,0)
        test Y.b:getorigin() == point.from(0,0,0)
        --Y
        test Y.a:vol()==1 and Y.b:vol()==1
        test Y.a:getorigin() == point.from(0,0,0)
        test Y.b:getorigin() == point.from(0,0,0)
    end

    testset "evaluation" do
        --X
        test X({0,0,0},{}) == X:getorigin()
        test X({1,1,1},{}) == point.from(1,1,1)
        --Y
        test Y({0,0,0},{}) == Y:getorigin()
        test Y({1,1,1},{}) == point.from(1,1,1)
    end
end






    --check if direction 'i' is a singular direction
    terra hypercube:setorigin(o : point)
        for i = 0, N do
            self.I[i]:setorigin(o(i))
        end
    end

    terra hypercube:getorigin()
        var o : point
        for i = 0, N do
            o(i) = self.I[i]:getorigin()
        end
        return o
    end

    terra hypercube:getcenter()
        var c : point
        for i = 0, N do
            c(i) = self.I[i].center
        end
        return c
    end

    --check if direction 'i' is a singular direction
    terra hypercube:issingulardir(i : int)
        return self.I[i]:vol()==0
    end

    --dimension is computed at runtime as the dimension of the embedding
    --minus the singular dimensions
    terra hypercube:dim()
        var d = N
        for i = 0, N do
            if self:issingulardir(i) then
                d = d - 1
            end
        end
        return d
    end

    --compute the length, area, volume, etc. 'vol' of a point
    --is defined as '1'
    terra hypercube:vol()
        var v = 1.0
        for i = 0, N do
            if not self:issingulardir(i) then
                v = v * self.I[i]:vol()
            end
        end
        return v
    end

    terra hypercube:barycentriccoord(p : point)
        var x : ntuple(T, N)
        var ptr = [&T](&x)
        var I : interval
        var k : uint8
        var D = self:dim()
        for i = 0, D do
            k = self.perm[i]
            I = self.I[k]
            ptr[i] = I:barycentriccoord(p(k))
        end
        for i = D, N do
            ptr[i] = 0
        end
        return x
    end

    local terra eval(self : &hypercube, x : &T)
        var d = self:dim()
        var p : point
        for i = 0, d do
            var k = self.perm[i]
            p(k) = self.I[k](x[i])
        end
        for i = d, N do
            var k = self.perm[i]
            p(k) = self.I[k].center
        end
        return p
    end

    hypercube.metamethods.__apply = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        local D = #args
        --try a 'point' or a 'tuple'
        if D==1 then
            local p = args[1]
            local typ = p.tree.type
            --try 'point'
            if typ==point then
                return `eval(&self, [&T](&[ p ]))
            else 
                --try n-tuple
                if typ.convertible == "tuple" then
                    local n = #typ.entries
                    local x = symbol(T[n])
                    local fill_array = terralib.newlist{}
                    for k = 0, n-1 do
                        fill_array:insert(quote [x][k] = p.["_"..tostring(k)] end)
                    end
                    return quote
                        err.assert(self:dim() == [n])
                        var [ x ]
                        [fill_array]
                    in
                        eval(&self, [&T](&x))
                    end
                end
            end
        end
        --try list of arguments
        local x = symbol(T[D])
        local fill_array = terralib.newlist{}
        for k,v in ipairs(args) do
            assert(v.tree.type:isprimitive() or "Not a primitive type")
            fill_array:insert(quote [x][k-1] = [ v ] end)
        end
        return quote
            var [ x ]
            err.assert(self:dim()==D)
            [ fill_array ]
        in
            eval(&self, [&T](&x))
        end
    end))

    --testing for equality of intervals
    hypercube.metamethods.__eq = terra(self : hypercube, other : hypercube) : bool
        for k = 0, N do
            if self.I[k] ~= other.I[k] then
                return false
            end
        end
        return true
    end

    hypercube.metamethods.__ne = terra(self : hypercube, other : hypercube) : bool
        return not (self==other)
    end
    
    hypercube.staticmethods.intersect = terra(self : hypercube, other : hypercube)
        var J : interval[N]
        for k = 0, N do
            var Z = interval.intersect(self.I[k], other.I[k])
            if Z:isempty() then
                return hypercube{J, false} --return empty (invalid) hypercube
            end
            J[k] = Z
        end
        return hypercube.new(J)
    end

    hypercube.metamethods.__mul = terra(a : hypercube, b : hypercube)
        var cube = hypercube.intersect(a, b)
        if not cube:isempty() then
            for i = 0, N do
                if a:issingulardir(i) and not b:issingulardir(i) then
                    cube.I[i] = b.I[i]
                elseif b:issingulardir(i) and not a:issingulardir(i) then
                    cube.I[i] = a.I[i]
                elseif a:issingulardir(i) and b:issingulardir(i) then
                    cube.I[i] = a.I[i]
                else
                    -- not a valid product - return with empty hypercube
                    cube.valid = false
                    return cube
                end
            end
        end
        cube:setperm()
        return cube
    end

    hypercube.metamethods.__div = terra(a : hypercube, b : hypercube)
        var cube = hypercube.intersect(a, b)
        var p = cube:getorigin()
        if cube == b then
            for i = 0, N do
                if b:issingulardir(i) then
                    cube.I[i] = a.I[i]
                else
                    cube.I[i].reach = 0
                    cube.I[i].center = cube.I[i].origin
                end
            end
        else
            cube.valid = false
        end
        cube:setperm()
        return cube
    end


    local intersection = terralib.overloadedfunction("intersection")
local T = double
for N = 2, 3 do
    for D1 = 0, N do
        for D2 = 0, N do
            local cube_d1, cube_d2 = Hypercube(T,N,D1), Hypercube(T,N,D1)
            intersection:adddefinition(
                terra(a : cube_d1, b : cube_d2)
                    var J : interval[N]
                    var ka : uint8, kb : uint8
                    var Z : interval
                    escape
                        for i = 0, N-1 do
                            emit quote
                                ka, kb = a.perm[i], b.perm[i] 
                                Z = intersect(a.I[ka], b.I[kb])
                                if Z:isempty() then
                                    return hypercube{J, false} --return empty (invalid) hypercube
                                end
                                J[k] = Z
                            end
                        end
                    end
                end
            )
        end
    end
end
