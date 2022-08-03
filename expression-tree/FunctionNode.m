classdef FunctionNode < ExpressionNode
   methods
       function infix = toInfix(obj)
            infix = obj.children{1}.toInfix();
            for i=2:numel(obj.children)
                infix = infix + "," + obj.children{i}.toInfix();
            end
            infix = string(obj.val) + "(" + infix + ")";
       end
   end
end