clear all
close all
clc

%%%
% This script is for loading and visualising two stream instability
% simulations in position/velocity phase space.
%%%

np = 8192; % total number of simulated electrons

times = 1:0.5:80; % list of times to be pulled into movie

figure('color','white')

% intialise video object
vidObj = VideoWriter('twoStream_Movie.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 5;
open(vidObj);

% loop through solution data and build movie
for t = 1:length(times)
    time = '../' + string(times(t));
    x = importPosition(time + '/lagrangian/pd/positions', 20,np+20);
    Ux = importVelocity(time + '/lagrangian/pd/U', 20,np+20);
    typeID = importTypeID(time + '/lagrangian/pd/typeId', 20,np+20);

    x_e1 = [];
    x_e2 = [];
    Ux_e1 = [];
    Ux_e2 = [];

    for i = 1:length(x)
        if(typeID(i) == 0)
            x_e1(end+1) = x(i);
            Ux_e1(end+1) = Ux(i);
        elseif(typeID(i) == 1)
            x_e2(end+1) = x(i);
            Ux_e2(end+1) = Ux(i);
        end
    end

    plot(x_e1',Ux_e1','.r')
    hold on
    plot(x_e2',Ux_e2','.b')
    hold off
    xlabel('x')
    ylabel('U_x')
    legend(['+U_x';'-U_x'])
    axis([0 6.221 -3 3])

    set(gca, 'nextplot','replacechildren');
    writeVideo(vidObj, getframe(gcf));
end

% close objects (good practise)
close(gcf)
close(vidObj);
