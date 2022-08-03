classdef BinaryNode < ExpressionNode
   methods
       function infix = toInfix(obj)
           infix = "(" + obj.children{1}.toInfix() + string(obj.val) + obj.children{2}.toInfix() + ")";
       end
   end
end