classdef VariableExpressionTree < ExpressionTree
    properties
        num_variables
    end
    methods (Access = protected)
        function [infix, counter] = variableToInfix(obj, node, counter)
            if isa(node, "FunctionNode")
                [infix, counter]= obj.variableToInfix(node.children{1}, counter);
                for i=2:numel(node.children)
                    [tmp_infix, counter]= obj.variableToInfix(node.children{i}, counter);
                    infix = infix + "," + tmp_infix;
                end
                infix = string(node.val) + "(" + infix + ")";
            elseif isa(node, "BinaryNode")
                [tmp_infix1, counter]= obj.variableToInfix(node.children{1}, counter);
                [tmp_infix2, counter]= obj.variableToInfix(node.children{2}, counter);
                infix = "(" + tmp_infix1 + string(node.val) + tmp_infix2 + ")";
            elseif isa(node, "LeafNode")
                if isa(node, "VariableLeafNode")
                    counter = counter + 1;
                    infix = sprintf("%s(%i)", node.val, counter);
                else
                    infix = string(node.val);
                    if isnumeric(node.val) && node.val < 0
                        infix = "(" + infix + ")";
                    end   
                end
            end
        end
        function [node, counter] = toExpressionNode(obj, node, variables, counter)
            if isa(node, "LeafNode")
                if isa(node, "VariableLeafNode")
                    counter = counter + 1;
                    node.val = variables(counter);
                end
            else
                for i=1:numel(node.children)
                    [node.children{i}, counter] = obj.toExpressionNode(node.children{i}, variables, counter);
                end
            end
        end
        function counter = countVariableLeafNode(obj, node, counter)
            if isa(node, "VariableLeafNode")
                counter = counter + 1;
            elseif ~isa(node, "LeafNode")
                for i=1:numel(node.children)
                    counter = obj.countVariableLeafNode(node.children{i}, counter);
                end
            end
        end
    end
    methods 
        function obj = VariableExpressionTree(root, num_variables)
            obj@ExpressionTree(root);
            if nargin == 2
                 obj.num_variables = num_variables;
            end
        end
        function func = toFunc(obj, variables)
            if nargin == 2
                infix = "@(x)" + obj.toInfix(variables);
            else
                infix = "@(x,c)" + obj.toInfix();
            end
            func = str2func(infix);
        end
        function infix = toInfix(obj, variables)
            if nargin == 2
                tree = obj.toExpressionTree(variables);
                infix = tree.toInfix();
            else
                [infix, ~] = obj.variableToInfix(obj.root, 0);
            end
        end
        function tree = toExpressionTree(obj, variables)
            node = obj.toExpressionNode(obj.root, variables, 0);
            tree = ExpressionTree(node);
        end
        
        function variables = optimize(obj, optimizer_name, x, objective_func, max_iter)
%             variables = rand(1, obj.num_variables);
            func = obj.toFunc();
            c_func = @(c) func(x, c);
            obj_func = @(c) objective_func(c_func(c));
            if obj.num_variables == 0
                variables = [];
            elseif optimizer_name == "particleswarm"
                options = optimoptions('particleswarm','MaxIterations', max_iter, 'Display', 'off', 'SwarmSize', obj.num_variables+2);
                variables = particleswarm(obj_func, obj.num_variables, [], [], options);
            elseif optimizer_name == "levenberg_marquadt"
                variables = levenberg_marquadt(obj_func, obj.num_variables, 1e-3, max_iter);
            elseif optimizer_name == "none"
                variables = (rand(obj.num_variables) - 0.5) * 20;
            else
                variables = eval(sprintf("%s(obj_func, %i)", optimizer_name, obj.num_variables));
            end
        end
        function tree = optimized(obj, optimizer_name, x, objective_func, max_iter)
            variables = obj.optimize(optimizer_name, x, objective_func, max_iter);
            if any(ismissing(variables))
                variables(ismissing(variables)) = 0;
            end
            tree = obj.toExpressionTree(variables);
        end
        function num_veriables = getNumVariables(obj)
            num_veriables = obj.countVariableLeafNode(obj.root, 0);
        end
    end
   methods(Static)
       function obj = randomTree(terminate_rate, max_depth, operand_rate, operands, operator_weights, operator_names, operator_argc)
            operator_weights = operator_weights ./ sum(operator_weights);
            operator_weights_interval = [cumsum(operator_weights), inf];
            [node, variable_count] = VariableExpressionTree.randomNode(terminate_rate, max_depth, operand_rate, operands, ...
                operator_weights_interval, 1, 0, operator_names, operator_argc);
            obj = VariableExpressionTree(node, variable_count);
       end
       function [obj , variable_count] = randomNode(terminate_rate, max_depth, operand_rate, operands, operator_weights_interval, depth, variable_count, ...
               operator_names, operator_argc)
            if rand() < terminate_rate || depth>=max_depth
                if rand() < operand_rate
                    idx = ceil(rand()*numel(operands));
                    obj = ExpressionNode.newNode(operands(idx));
                else
                    variable_count  = variable_count + 1;
                    obj = VariableLeafNode("c");
                end
            else
                idx = find(rand() < operator_weights_interval, true, 'first');
                op = operator_names(idx);
                argc = operator_argc(idx);
                children = cell(1, argc);
                is_leaf = zeros(1, argc);
                for i=1:argc
                    [children{i}, variable_count] = VariableExpressionTree.randomNode(terminate_rate, max_depth, operand_rate, operands, ...
                            operator_weights_interval, depth+1, variable_count, operator_names, operator_argc);
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
       
       
   end
end 