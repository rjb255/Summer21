clear
comp = {};
stage = "1";
comp{1, 1} = readmatrix("epsilon_linspace.csv");
comp{2, 1} = readmatrix("epsilon_fmin.csv");
comp{3, 1} = readmatrix("epsilon_aLearning.csv");
x = repmat(2:size(comp{1,1}, 1)+1, size(comp{1,1}, 2), 1)';
xp = (2:size(comp{1,1}, 1)+1)';
colour = ["k", "b", "r"];
marks = [".", "x", "+"];

for i = 1:size(comp, 1)
    comp{i, 2} = mean(comp{i,1}, 2);
    comp{i, 3} = sqrt(var(comp{i,1}, 0, 2));
    comp{i, 4} = abs(mean(comp{i,1}, 2));
end
a = figure;%("visible", "off");

t = tiledlayout(2,2, "TileSpacing", 'compact');

nexttile(3, [1, 2])
plot(xp, zeros(size(xp)), 'k:', "HandleVisibility", "off")
hold on
for i = 1:size(comp, 1)
    plot(xp, comp{i, 1}, colour(i) + marks(i), "HandleVisibility", "off")
    plot(xp, comp{i, 2}, colour(i), "HandleVisibility", "on")
end
hold off
update()
leg = legend("Greatest Uncertainty", "fminbound", "Problem Specific", "Location", "southoutside", "Orientation", "horizontal");
set(leg, 'box', 'off', 'Interpreter', 'latex')
ylabel("$\varepsilon$", "Interpreter", "latex")

nexttile(1)
for i = 1:size(comp, 1)
    plot(xp, comp{i, 3}, colour(i), "HandleVisibility", "on")
    hold on
end
hold off
update()
ylabel("$\sigma_\varepsilon$", "Interpreter", "latex")

nexttile
for i = 1:size(comp, 1)
    plot(xp, comp{i, 4}, colour(i), "HandleVisibility", "on")
    hold on
end
hold off
update()
ylabel("$\left|\overline{\varepsilon}\right|$", "Interpreter", "latex")

a.Position = [10 10 900 900];
exportgraphics(t, "../LaTeXLearning\Version1\phd-thesis-template-2.4\Chapter2\Figs\Vector\comparison"+stage+".pdf", 'ContentType', 'vector')
exportgraphics(t, "../LaTeXLearning\Version1\phd-thesis-template-2.4\Chapter2\Figs\PDF\comparison"+stage+".pdf", 'ContentType', 'vector')

function update()
    xlabel("Number of Samples", "Interpreter", "latex")
    xlim([0,25])
    ax = gca;
    set(ax, "TickLabelInterpreter", "latex")
end