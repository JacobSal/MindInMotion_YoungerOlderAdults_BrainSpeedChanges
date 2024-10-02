% This code adapted from std_topoplot
function [STUDY] = std_topoplot_CL(STUDY,cls,mode)

icadefs;
if isempty(cls)
    cls = 3:length(STUDY.cluster); % plot all clusters in STUDY    
end
if strcmpi(mode,'apart')
    for clus = 1: length(cls) % For each cluster requested
        len = length(STUDY.cluster(cls(clus)).comps);
        if len > 0 % A non-empty cluster 
            h_topo = figure;
            rowcols(2) = ceil(sqrt(len + 4)); rowcols(1) = ceil((len+4)/rowcols(2));
            clusscalp = STUDY.cluster(cls(clus));
            ave_grid = clusscalp.topo;
            tmp_ave = ave_grid;
            tmp_ave(find(isnan(tmp_ave))) = 0; % remove NaN values from grid for later correlation calculation.  
            for k = 1:len
                abset   = STUDY.datasetinfo(STUDY.cluster(cls(clus)).sets(1,k)).index;
                subject = STUDY.datasetinfo(STUDY.cluster(cls(clus)).sets(1,k)).subject;
                comp = STUDY.cluster(cls(clus)).comps(k);
                [Xi,Yi] = meshgrid(clusscalp.topoy,clusscalp.topox);                     
                scalpmap = squeeze(clusscalp.topoall{k}); % already correct polarity
                if k <= rowcols(2) - 2 %first sbplot row
                    figure(h_topo);
                    sbplot(rowcols(1),rowcols(2),k+2) , 
                    toporeplot(scalpmap, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off','xsurface', Xi, 'ysurface', Yi );
                    title([subject '/' 'IC' num2str(comp)   ], 'interpreter', 'none');
                    colormap(DEFAULT_COLORMAP);
                    %waitbar(k/(len+1),h_wait)
                else %other sbplot rows
                    figure(h_topo)
                    sbplot(rowcols(1),rowcols(2),k+4) , 
                    toporeplot(scalpmap, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off','xsurface', Xi, 'ysurface', Yi );
                    title([subject '/' 'IC' num2str(comp)], 'interpreter', 'none');
                    colormap(DEFAULT_COLORMAP);
                    %waitbar(k/(len+1),h_wait)
                end
            end
            figure(h_topo)
            sbplot(rowcols(1),rowcols(2),[1 rowcols(2)+2 ]) ,
            toporeplot(ave_grid, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off');
            title([ STUDY.cluster(cls(clus)).name ' (' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(clus)).comps)),' ICs)']);
            %title([ STUDY.cluster(cls(clus)).name ' average scalp map, ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) 'Ss']);
            set(gcf,'Color', BACKCOLOR);
            colormap(DEFAULT_COLORMAP);
            %waitbar(1,h_wait)
            %delete(h_wait)
            orient tall  % fill the figure page for printing
            axcopy
        end % Finished one cluster plot 
    end   % Finished plotting all clusters
end

% Plot clusters centroid maps
if strcmpi(mode, 'together') 
    len = length(cls);
    rowcols(2) = ceil(sqrt(len)); rowcols(1) = ceil((len)/rowcols(2));
%     figure
    for k = 1:len 
       if len ~= 1
            sbplot(rowcols(1),rowcols(2),k)  
       end
        tmpcmap = colormap(DEFAULT_COLORMAP);
        toporeplot(STUDY.cluster(cls(k)).topo, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off','colormap', tmpcmap);
        title([ STUDY.cluster(cls(k)).name ' (' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(k)).comps)),' ICs)' ]);
        colormap(DEFAULT_COLORMAP);
        %title([ STUDY.cluster(cls(k)).name ', ' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ]);
    end
    if len ~= 1
        maintitle = 'Average scalp map for all clusters';
        a = textsc(maintitle, 'title'); 
        set(a, 'fontweight', 'bold');
        set(gcf,'name', maintitle);
    else
        title([ STUDY.cluster(cls(k)).name ' (' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(k)).comps)),' ICs)']);
        set(gcf,'name',['Scalp map of ' STUDY.cluster(cls(k)).name ' (' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) ' Ss, ' num2str(length(STUDY.cluster(cls(k)).comps)),' ICs)']);
        %title([ STUDY.cluster(cls(k)).name ' scalp map, ' num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ]);
    end
    set(gcf,'Color', BACKCOLOR);

    orient tall
    axcopy
end % Finished 'apart' plotting mode