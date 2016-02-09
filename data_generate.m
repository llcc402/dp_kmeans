function data = data_generate()

theta = rand(1000, 1) * 2 * pi;
r = 1;
epsilon = 0.2;

x = r * cos(theta) + epsilon * rand(1000, 1);
y = r * sin(theta) + epsilon * rand(1000, 1);

data1 = [x, y];

r = 2;
theta = rand(1000, 1) * 2 * pi;
x = r * cos(theta) + epsilon * rand(1000, 1);
y = r * sin(theta) + epsilon * rand(1000, 1);

data2 = [x, y];

data = [data1; data2];

plot(data(:,1), data(:, 2), '.')

end