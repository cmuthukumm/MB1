function Conclude1(DisplayFlag, OPTIONS, Population, nLegal, MinCost, AvgCost)



if DisplayFlag
    % Count the number of duplicates
    NumDups = 0;
    for i = 1 : OPTIONS.popsize
        Chrom1 = sort(Population(i).chrom);
        for j = i+1 : OPTIONS.popsize
            Chrom2 = sort(Population(j).chrom);
            if isequal(Chrom1, Chrom2)
                NumDups = NumDups + 1;
            end
        end
    end  
    disp([num2str(NumDups), ' duplicates in final population.']);
    disp([num2str(nLegal), ' legal individuals in final population.']);
    Chrom = sort(Population(1).chrom);
    disp(['Best chromosome = ', num2str(Chrom)]); 
   figure(2);
    plot((0:OPTIONS.Maxgen), MinCost, '-*r');
    hold on
    grid on;
    plot((0:OPTIONS.Maxgen), AvgCost, '-ob');
    xlabel('Iterations');
    ylabel('Fitness');
    legend('BO','MBO');
    title('Convergence Plot');
    hold off
end
% return;