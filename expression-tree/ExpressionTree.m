classdef ExpressionTree
    properties (Constant)
        binary_operator_names = ["+", "-", ".*", "./", ".^"];
        function_operator_names = [ "sin", "cos"];
        operator_names = cat(2, ExpressionTree.binary_operator_names, ExpressionTree.function_operator_names);
        
        operators = containers.Map( ...
            ExpressionTree.operator_names, ...
            {@(a,b)a+b, @(a,b)a-b, @(a,b)a*b, @(a,b)a/b, @(a,b)a^b, @(a)sin(a), @(a)cos(a)} ...
        );
        
        priority = containers.Map( ...
            ['+', '-', '*', '/', '^'], ...
            [ 1 ,  1 ,  2 ,  2 ,  3 ] ...
        );

        operator_argc = containers.Map( ...
            ExpressionTree.operator_names, ...
            [2, 2, 2, 2, 2, 1, 1] ...
        );
    end
    properties
        root
    end
    methods 
        function obj = ExpressionTree(root)
            if nargin == 1
                 obj.root = root;
            else
                 obj.root = 0;
            end
        end
        function func = toFunc(obj)
            func = str2func("@(x)" + obj.toInfix());
        end
        function infix = toInfix(obj)
           infix = obj.root.toInfix();
        end
        function txt = toLatex(obj, num_classes)
            infix = obj.toInfix();
            args = [];
            for i=1:num_classes
                infix = strrep(infix, sprintf("x(:,%i)", i), sprintf("x_%i", i));
                args = [args, sprintf("x_%i", i)];
            end
            args = strjoin(args, ",");
            func = str2func(sprintf("@(%s)%s", args, infix));
            func = sym(func);
            txt = latex(func);
        end
        function depth = getDepth(obj)
            depth = obj.root.getDepth();
        end
        function count = getNodeCount(obj)
            count = obj.root.getNodeCount();
        end
        function ptrs = getPtrAtDepth(obj, node, depth)
           depth = depth-1;
           if depth == 1
               ptrs = num2cell(1:numel(node.children));
           elseif depth >= 2 && ~isa(node, "LeafNode")
               ptrs = cell(1, numel(node.children));
               for i=1:numel(node.children)
                   ptr = obj.getPtrAtDepth(node.children{i}, depth);
                   ptr = cellfun(@(n) [i, n], ptr, 'UniformOutput', false);
                   ptrs{i} = ptr;
               end
               ptrs = cat(2, ptrs{:});
           elseif depth == 0
               ptrs = {[]};
           else
               ptrs = {};
           end
        end
        function node = getNodeAt(obj, node, ptr)
            if numel(ptr)>0 && ptr(end) == 0
                node = NullNode(0);
            else
                for p=ptr
                    node = node.children{p};
                end
            end
        end
        function node = setNodeAt(obj, node, ptr, newnode)
            if numel(ptr) == 1
                node.children{ptr(1)} = newnode;
            elseif numel(ptr) == 0
                node = newnode;
            elseif isa(newnode, "NullNode") 
%                 node = node;
            else
                node.children{ptr(1)} = obj.setNodeAt(node.children{ptr(1)}, ptr(2:end), newnode);
            end
        end
   end
   methods(Static)
       function obj = randomTree(terminate_rate, max_depth, operand_rate, operands, operator_weights, constant_range)
            obj = ExpressionTree();
            operator_weights = operator_weights ./ sum(operator_weights);
            operator_weights_interval = [cumsum(operator_weights), inf];
            node = ExpressionTree.randomNode(terminate_rate, max_depth, operand_rate, operands, operand_weights_interval, constant_range, 1);
            obj.root = node;
       end
       function obj = randomNode(terminate_rate, max_depth, operand_rate, operands, operator_weights_interval, constant_range, depth)
            if rand() < terminate_rate || depth>=max_depth
                if rand() < operand_rate
                    idx = ceil(rand()*numel(operands));
                    obj = ExpressionNode.newNode(operands(idx));
                else
                    obj = ExpressionNode.newNode(rand()*(constant_range(2) - constant_range(1)) + constant_range(1));
                end
            else
                idx = find(rand() < operator_weights_interval, true, 'first');
                op = ExpressionTree.operator_names(idx);
                argc = ExpressionTree.operator_argc(op);
                children = cell(1, argc);
                is_leaf = zeros(1, argc);
                for i=1:argc
                    children{i} = ExpressionTree.randomNode(terminate_rate, max_depth, operand_rate, operands, operator_weights_interval, constant_range, depth+1);
                    is_leaf(i) = isa(children{i}, "LeafNode");
                end
                % ensure there are at least one operand in formula
                if all(is_leaf)
                    is_const = zeros(1, argc);
                    for i=1:argc
                        is_const(i) = isnumeric(children{i}.val) || ~any(ismember(operands, children{i}.val));
                    end
                    if all(is_const)
                        idx = ceil(rand()*argc);
                        children{idx} = ExpressionNode.newNode(operands(ceil(rand()*numel(operands))));
                    end
                end
                obj = ExpressionNode.newNode(op, children);
            end
       end
       
       function obj  = fromInfix(infix, variables)
           if isa(infix, 'string')
               infix = infix.char;
           end
           assert(isa(infix, 'char'));
           if nargin == 1
               variables = [];
           end

%            sInfix = convertCharsToStrings(infix);
           flag = 0;
           operator_keys = ExpressionTree.priority.keys();
           tmp_str = "";
           for i=1:numel(infix)
               ch = infix(i);
               if regexp(ch,'[a-zA-Z]') == 1
                   tmp_str = tmp_str + ch;
                   flag = 1;
               elseif ch == '('
                   if flag == 1 % sin(
                       
                   end
              
               elseif ch == ')'
                   
               elseif contains(operator_keys{1}, ch)
               elseif any(ismember(ExpressionTree.func_name, tmp_str))
               end
           end
       end
   end
end 