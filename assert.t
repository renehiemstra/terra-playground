local C = terralib.includecstring [[
    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    #include <stdarg.h>
]]

local S = {}

S.assert = macro(function(condition, message)
    local loc = condition.tree.filename..":"..condition.tree.linenumber
    local print_message = quote end
    if message then
        print_message = quote C.printf("%s\n", message) end
    end
    return quote
	    if not condition then
            C.printf("%s: assertion failed!\n", loc)
            [print_message]
    	    C.abort()
	    end -- if
    end -- quote
end) -- macro

S.error = macro(function(expr)
    local loc = expr.tree.filename..":"..expr.tree.linenumber
    return quote
	    C.printf("%s: %s\n", loc, expr)
	    C.abort()
    end -- quote
end) -- macro

return S

