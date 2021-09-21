clear
comp = {};
comp{1, 1} = readmatrix("epsilon_linspace.csv");
comp{2, 1} = readmatrix("epsilon_fmin.csv");
comp{3, 1} = readmatrix("epsilon_aLearning.csv");

colour = ["k", "b", "r"];
marks = [".", "x", "+"];

for i = 1:size(comp, 1)
    comp{i, 2} = mean(comp{i,1}, 2);
    comp{i, 3} = sqrt(var(comp{i,1}, 0, 2));
    comp{i, 4} = sqrt(mean(comp{i,1}.^2, 2));
end

subplot(2,1,2)
plot([0, 15], [0, 0], 'k:', "HandleVisibility", "off")
hold on
for i = 1:size(comp, 1)
    plot(comp{i, 1}, colour(i) + marks(i), "HandleVisibility", "off")
    plot(comp{i, 2}, colour(i), "HandleVisibility", "on")
end
hold off
update()
ylabel("$\varepsilon$", "Interpreter", "latex")

subplot(2,2,1)
for i = 1:size(comp, 1)
    plot(comp{i, 3}, colour(i), "HandleVisibility", "on")
    hold on
end
hold off
update()
ylabel("$\sigma_\varepsilon$", "Interpreter", "latex")

subplot(2,2,2)
for i = 1:size(comp, 1)
    plot(comp{i, 4}, colour(i), "HandleVisibility", "on")
    hold on
end
hold off
update()
ylabel("$\sqrt{\overline{\varepsilon^{2}}}$", "Interpreter", "latex")

exportgraphics(gcf, "../LaTeXLearning\Version1\phd-thesis-template-2.4\Chapter2\Figs\Vector\comparison1.pdf", 'ContentType', 'vector')
exportgraphics(gcf, "../LaTeXLearning\Version1\phd-thesis-template-2.4\Chapter2\Figs\PDF\comparison1.pdf", 'ContentType', 'vector')

function update()
    xlabel("$x$", "Interpreter", "latex")
    legend("Linear", "fminbound", "Active Learning")
end

