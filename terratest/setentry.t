local C = terralib.includec("stdio.h")
import "terratest"

struct A{
    b : int
}

testenv "setentry" do

    testset "setentry" do
        terracode
            var a = A{1}
            a.b = 2
        end
        test a.b == 2
    end

end