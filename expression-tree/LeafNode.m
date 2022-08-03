classdef LeafNode < ExpressionNode
    methods
        function infix = toInfix(obj)
            infix = string(obj.val);
            if isnumeric(obj.val) && obj.val < 0
                infix = "(" + infix + ")";
            end
        end
        function depth = getDepth(obj)
            depth = 1;
        end
        function count = getNodeCount(obj)
            count = 1;
        end
       
   end
end