TestTable = DataOut.AllMice.ClusterData%Up %Down %(~(ismember(DataOut.AllMice.ClusterData.Area, vecROI)),:);
Areas = unique(TestTable.Area);
Frequencies = [];
for i = 1:length(Areas)
    Frequencies = [Frequencies; sum(TestTable.Area == Areas(i))];
end
Areas = categorical(Areas);

figure;
bar(Areas, Frequencies);
% set(gca,'xticklabel',Areas);
ylim([0 200]);

