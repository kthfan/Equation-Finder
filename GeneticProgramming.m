classdef GeneticProgramming
   properties
      population_size
      selection_rate
      individual_pool
      optimized_individual_pool
      num_classes
      fitness_func
      mutation_rate
      history
      terminate_rate
      max_depth
      operands
      operand_rate
      constant_range
      operator_weights
      best_equation
      best_fitness
      constant_optimizer
      operator_names
      operator_argc
      constant_optimize_iter
      last_execution_time
   end
   methods
       function obj = GeneticProgramming(varargin)
           obj.last_execution_time = 0;
           obj.population_size = 100;
           obj.selection_rate = 0.5;
           obj.best_fitness = inf;
           obj.individual_pool = cell(0);
           obj.optimized_individual_pool = cell(0);
           obj.terminate_rate = 0.05;
           obj.max_depth = 5;
           obj.operand_rate = 0.75;
           obj.operator_weights = [5, 5, 5, 5, 3, 1, 1];
           obj.constant_range = [-10, 10];
           obj.fitness_func = @(pred_y, true_y) sum(abs(pred_y - true_y).^2);
           obj.history = containers.Map(["fitness", "equation"], {[], {}});
           obj.mutation_rate = 0.1;
           obj.constant_optimizer = "levenberg_marquadt";
           obj.operator_names = ExpressionTree.operator_names;
           obj.operator_argc = [2, 2, 2, 2, 2, 1, 1];
           obj.constant_optimize_iter = 5;
           i = 1;
           while i <= numel(varargin)
               switch varargin{i}
                   case 'PopulationSize'
                       obj.population_size = varargin{i+1};
                   case 'SelectionRate'
                       obj.selection_rate = varargin{i+1};
                   case 'TreeTerminateRate'
                       obj.terminate_rate = varargin{i+1};
                   case 'TreeMaxDepth'
                       obj.max_depth = varargin{i+1};
                   case 'TreeOperandRate'
                       obj.operand_rate = varargin{i+1};
                   case 'ConstantRange'
                       obj.constant_range = varargin{i+1};
                   case 'Classes'
                       obj.num_classes = varargin{i+1};
                   case 'FitnessFunction'
                       obj.fitness_func = varargin{i+1};
                   case 'OperatorWeights'
                       obj.operator_weights = varargin{i+1};
                   case 'MutationRate'
                       obj.mutation_rate = varargin{i+1};
                   case 'ConstantOptimizer'
                       obj.constant_optimizer = varargin{i+1};
                   case 'Operators'
                       obj.operator_names = varargin{i+1};
                   case 'OperatorParameters'
                       obj.operator_argc = varargin{i+1};
                   case 'ConstantOptimizeIterations'
                       obj.constant_optimize_iter = varargin{i+1};
               end
               i = i + 2;
           end
           
       end
       function obj = initialize(obj)
           % generate operands: x(:, 1), x(:, 2), ...
           obj.operands = arrayfun(@(i) sprintf("x(:,%i)", i), 1:obj.num_classes);
           
           % initialize ExpressionTree
           individual_pool = cell(1, obj.population_size);
           for i=1:obj.population_size
%                individual_pool{i} = ExpressionTree.randomTree(obj.terminate_rate, obj.max_depth, obj.operand_rate, obj.operands, operand_weights, obj.constant_range);
               individual_pool{i} = VariableExpressionTree.randomTree(obj.terminate_rate, obj.max_depth, ...
                   obj.operand_rate, obj.operands, obj.operator_weights, obj.operator_names, obj.operator_argc);
           end
           obj.individual_pool = individual_pool;
       end
       function [node1, offset] = getNodeAtOffset(obj, node, offset)
           offset = offset - 1;
           node1 = node;
           if offset ~= 0 && ~isa(node, "LeafNode")
               for i=1:numel(node.children)
                   [node1, offset] = obj.getNodeAtOffset(node.children{i},offset);
                   if offset==0
                       break;
                   end
               end
           end
       end
       
       function nodes = getNodesAtDepth(obj, node, depth)
           depth = depth-1;
           if depth == 1
               nodes = node.children;
           elseif depth >= 2 && ~isa(node, "LeafNode")
               nodes = cellfun(@(n) obj.getNodesAtDepth(n, depth-1), node.children, 'UniformOutput', false);
               nodes = cat(2, nodes{:});
           elseif depth == 0 && isa(node, "LeafNode")
               nodes = {node};
           else
               nodes = {};
           end
       end
       function [x, y] = getCutDepth(obj, a, b, M)
           lower = max(1-b, -M+a);
           upper = min(-1+a, M-b);
           r = round(rand()*(upper-lower)+lower);
           la = max(1, r+1);
%            lb = max(1-r, 1);
           ua = min(a, b+r);
%            ub = min(b, a-r);
           x = round(rand()*(ua-la)+la);
%          y = round(rand()*(ub-lb)+lb);
           y = x-r;
           assert(x-y==r)
       end 
       function [tree1, tree2] = crossoverTree(obj, tree1, tree2, max_depth)
%            node_count1 = tree1.getNodeCount();
           depth1 = tree1.getDepth();
           depth2 = tree2.getDepth();
           % generate depths to crossover
           [cut_depth1, cut_depth2] = obj.getCutDepth(depth1, depth2, max_depth);
           cut_depth1 = depth1 - cut_depth1 + 1; cut_depth2 = depth2 - cut_depth2 + 1;
           cut_ptrs1 = tree1.getPtrAtDepth(tree1.root, cut_depth1);
           cut_ptrs2 = tree2.getPtrAtDepth(tree2.root, cut_depth2);

           cut_index1 = ceil(numel(cut_ptrs1)*rand()); % index of childern nodes that contain cut_nodes1
           cut_index2 = ceil(numel(cut_ptrs2)*rand()); % index of childern nodes that contain cut_nodes2

           % check whether to switch
           switch1 = false; switch2 = false;
           if cut_index1 ~= 0
               cut_ptr1 = cut_ptrs1{cut_index1};
               cut_node1 = tree1.getNodeAt(tree1.root, cut_ptr1);
               switch1 = true;
           end
           if cut_index2 ~= 0
               cut_ptr2 = cut_ptrs2{cut_index2};
               cut_node2 = tree2.getNodeAt(tree2.root, cut_ptr2);
               switch2 = true;
           end
           if switch1
               tree2.root = tree2.setNodeAt(tree2.root, cut_ptr2, cut_node1);
           end
           if switch2
               tree1.root = tree1.setNodeAt(tree1.root, cut_ptr1, cut_node2);
           end
           tree1.num_variables = tree1.getNumVariables();
           tree2.num_variables = tree2.getNumVariables();
       end
       function [obj, selected_individual] = tournament_selection(obj, x, y)
           selection_size = round(obj.population_size * obj.selection_rate);
           func_list = cell(1, obj.population_size);
           fitness_vals = zeros(1, obj.population_size);
           % calculate fitness of equations
           for i=1:obj.population_size
               func_list{i} = obj.optimized_individual_pool{i}.toFunc();
               fitness_vals(i) = obj.fitness_func(func_list{i}(x), y);
           end
           
           % sort individuals by fitness
           [fitness_vals, order] = sort(fitness_vals);
           obj.individual_pool = obj.individual_pool(order);
           obj.optimized_individual_pool = obj.optimized_individual_pool(order);
           selected_individual = obj.individual_pool(1:selection_size);

           % record history
           if fitness_vals(1) <= obj.best_fitness
               obj.best_equation = func_list{order(1)};
               obj.best_fitness = fitness_vals(1);
           end
           
           obj.history("fitness") = [obj.history("fitness"), obj.best_fitness];
           history_equations = obj.history("equation");
           history_equations{end+1} = obj.best_equation;
           obj.history("equation") = history_equations;
       end
       function descendants = crossover(obj, selected_individual)
           mutation_size = round(obj.population_size * obj.mutation_rate);
           descendant_size = obj.population_size - round(obj.population_size * obj.selection_rate) - mutation_size;
           % if selected_individual is odd, transform into even
           if mod(numel(selected_individual), 2) == 1
               selected_individual{end + 1} = selected_individual{1};
           end
           selection_size = numel(selected_individual);
           
           selected_individual = selected_individual(randperm(selection_size)); % shuffle
           
           if descendant_size <= selection_size
               tree_pairs1 = cell(1, ceil(descendant_size/2));
               tree_pairs2 = cell(1, ceil(descendant_size/2));
               for i=1:ceil(descendant_size/2)
                   [tree_pairs1{i}, tree_pairs2{i}] = obj.crossoverTree(selected_individual{2*i-1}, selected_individual{2*i}, obj.max_depth);
               end
               descendants = cat(2, tree_pairs1, tree_pairs2);
               descendants = descendants(1:descendant_size);
           else
               z_iter = floor(descendant_size / selection_size);
               z_remainder = mod(descendant_size, selection_size);
               if z_remainder==0
                   descendants = cell(1, z_iter);
               else
                   descendants = cell(1, z_iter+1);
               end
               for z=1:z_iter
                   tree_pairs1 = cell(1, ceil(selection_size/2));
                   tree_pairs2 = cell(1, ceil(selection_size/2));
                   for i=1:ceil(selection_size/2)
                       [tree_pairs1{i}, tree_pairs2{i}] = obj.crossoverTree(selected_individual{2*i-1}, selected_individual{2*i}, obj.max_depth);
                   end
                   z_descendants = cat(2, tree_pairs1, tree_pairs2);
                   z_descendants = z_descendants(1:selection_size);
                   descendants{z} = z_descendants;
               end
               if z_remainder~=0
                   tree_pairs1 = cell(1, ceil(z_remainder/2));
                   tree_pairs2 = cell(1, ceil(z_remainder/2));
                   for i=1:ceil(z_remainder/2)
                       [tree_pairs1{i}, tree_pairs2{i}] = obj.crossoverTree(selected_individual{2*i-1}, selected_individual{2*i}, obj.max_depth);
                   end
                   z_descendants = cat(2, tree_pairs1, tree_pairs2);
                   descendants{z_iter+1} = z_descendants(1:z_remainder);
               end
               descendants = horzcat(descendants{:});
           end 
       end
       function mutations = mutation(obj, selected_individual)
           mutation_size = round(obj.population_size * obj.mutation_rate);
           % if selected_individual is odd, transform into even
           if mod(numel(selected_individual), 2) == 1
               selected_individual{end + 1} = selected_individual{1};
           end
           selection_size = numel(selected_individual);

           z_iter = floor(mutation_size / (2*selection_size));
           z_remainder = mod(mutation_size, 2*selection_size);
           if z_remainder==0
               mutations = cell(1, z_iter);
           else
               mutations = cell(1, z_iter+1);
           end
           for z=1:z_iter
               tree_pairs1 = cell(1, selection_size);
               tree_pairs2 = cell(1, selection_size);
               for i=1:ceil(selection_size)
                   new_tree = VariableExpressionTree.randomTree(obj.terminate_rate, obj.max_depth, ...
                                obj.operand_rate, obj.operands, obj.operator_weights, obj.operator_names, obj.operator_argc);
                   [tree_pairs1{i}, tree_pairs2{i}] = obj.crossoverTree(selected_individual{i}, new_tree, obj.max_depth);
               end
               z_mutations = cat(2, tree_pairs1, tree_pairs2);
               mutations{z} = z_mutations;
           end
           if z_remainder~=0
               tree_pairs1 = cell(1, ceil(z_remainder/2));
               tree_pairs2 = cell(1, ceil(z_remainder/2));
               for i=1:ceil(z_remainder/2)
                   new_tree = VariableExpressionTree.randomTree(obj.terminate_rate, obj.max_depth, ...
                                obj.operand_rate, obj.operands, obj.operator_weights, obj.operator_names, obj.operator_argc);
                   [tree_pairs1{i}, tree_pairs2{i}] = obj.crossoverTree(selected_individual{i}, new_tree, obj.max_depth);
               end
               z_mutations = cat(2, tree_pairs1, tree_pairs2);
               mutations{z_iter+1} = z_mutations;
           end
           mutations = horzcat(mutations{:});
           mutations = mutations(1:mutation_size);

       end
       function obj = opimizeTrees(obj, x, y)
           objective_func = @(y_pred) obj.fitness_func(y, y_pred);
           optimized_individual_pool = cell(1, numel(obj.individual_pool));
           for i=1:numel(obj.individual_pool)
               optimized_individual_pool{i} = obj.individual_pool{i}.optimized(obj.constant_optimizer, x, objective_func, obj.constant_optimize_iter);
           end
           obj.optimized_individual_pool = optimized_individual_pool;
       end
       function obj = single_generation(obj, x, y, iter)
           start_time = now;
           obj = obj.opimizeTrees(x, y);
           [obj, selected_individual] = obj.tournament_selection(x, y);
           descendants = obj.crossover(selected_individual);
           mutations = obj.mutation(selected_individual);
           obj.individual_pool = cat(2, selected_individual, descendants, mutations);
           obj.last_execution_time = (now - start_time) * 1e5;
       end
       function obj = fit(obj, x, y, max_generations)
           for iter=1:max_generations
               obj = obj.single_generation(x, y, iter);
           end
       end
       function y = predict(obj, x)
           y = obj.best_equation(x);
       end
   end
end