classdef NullNode < ExpressionNode
    methods
        function infix = toInfix(obj)
            infix = "";
        end
        function depth = getDepth(obj)
            depth = 0;
        end
        function count = getNodeCount(obj)
            count = 0;
        end
       
   end
end