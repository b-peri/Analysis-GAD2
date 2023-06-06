for i = 1:16
    subplot(4,4,i);
    plot(sAP.sCluster(i).Waveform);
    title(["Cluster: "+i])
end