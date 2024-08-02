
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