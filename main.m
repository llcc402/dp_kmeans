data = data_generate();
tic;
z = dp_kmeans(data);
toc

m = max(z);
figure(2)
hold on
xlim([-3,3])
ylim([-3,3])
for i = 1:m
    try
        plot(data(z==i,1),data(z==i,2),'.')
    end
end