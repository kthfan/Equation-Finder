

addpath("./expression-tree");
addpath("./utils/optimizers");


GP = GeneticProgramming('Classes', 1, 'PopulationSize', 50, "TreeMaxDepth", 3, ...
            'SelectionRate', 0.213);
GP = GP.initialize();

x = linspace(0.001, 3, 50)';
y = cos(2*x);

GP = GP.fit(x, y, 20);

func = GP.history("equation");
func = func{1};

figure(1)
plot(x, func(x))
hold on
plot(x, y)
hold off







