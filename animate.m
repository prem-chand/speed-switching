%vars = dlmread('BiMASC_config.dta','\t');
disp('Animation corresponds to data in the workspace')

t       = simout.t;
pstk    = simout.stk;
ph      = simout.h;
pswk    = simout.swk;
pswf    = simout.swf;
phead   = simout.head;

vars = [t pstk ph pswk pswf phead];



torig = vars(:,1);
p     = vars(:,2:11);


% if 0
%   torig = logsout.FullState.time';
%   q = logsout.FullState.Data(:,1:7)';
%   pH = logsout.FullState.Data(:,8:9)';
%   foot_switch = logsout.FootSwitches.Data';
% end
% 
% if 1
%   torig = simout.t;
%   q = simout.x(1:7,:);
%   pH = simout.terminal_output.pH;
%   foot_switch = [ones(1,length(simout.t));zeros(1,length(simout.t))];
% end

RATE = 100;

%==========================================
% initialize the animation figure and axes
%==========================================

fig1 = figure;

set(0,'Units','pixels')
scnsize = get(0,'ScreenSize');

screen_width = scnsize(3);
screen_height = scnsize(4);

figure_x_limits = [-1.5 4];
figure_y_limits = [-.2 2.8];

% find the minimum scaling factor

figure_x_size = figure_x_limits(2) - figure_x_limits(1);
figure_y_size = figure_y_limits(2) - figure_y_limits(1);

xfactor = screen_width/figure_x_size;
yfactor = screen_height/figure_y_size;

if (xfactor < yfactor)
  screen_factor = 0.9*xfactor;
else
  screen_factor = 0.9*yfactor;
end

% calculate screen offsets
screen_x_offset = (screen_width - screen_factor*figure_x_size)/2;
screen_y_offset = (screen_height - screen_factor*figure_y_size)/2;

% draw figure and axes
set(fig1,'Position', [screen_x_offset screen_y_offset screen_factor*figure_x_size screen_factor*figure_y_size]);
set(fig1,'MenuBar', 'none');
axes1 = axes;
set(axes1,'XLim',figure_x_limits,'YLim',figure_y_limits);
set(axes1,'Position',[0 0 1 1]);
set(axes1,'Color','w');
set(axes1,'TickDir','out');
%sp1=size(p)
%[t, p] = even_sample(torig, p, RATE)
%sp2=size(p)
%str_lng(:,1)
%hhh=ceil(str_lng(:,1)*sp2/sp1)

%[t, pH] = even_sample(torig', pH', RATE);
%[t, foot_switch] = even_sample(torig', foot_switch', RATE);

p=p';
%pH=pH';
%foot_switch=foot_switch';
t=t';

box on;

MAKE_MOVIE = 0;
shift=[0 0];
%drawing the terrain


for i=1:1:size(p,2)
  %drawone(q(:,i),pH(:,i),foot_switch(:,i),axes1);
%   for j=size(str_lng,1)-1:-1:1
%       if i>str_lng(j,1)
%          shift=str_lng(j,2:3);
%          break;
%       end
%   end
  drawone(p(:,i),axes1,shift);
  s = sprintf('Hopping Fixed point (fwd speed in flight ~1.6 m/s)',t(i));
  text(-1.3,2.4,s,'FontAngle','italic','FontWeight','bold');
  drawnow;
  %if MAKE_MOVIE, M(i) = getframe; end
end

if MAKE_MOVIE
  movie2avi(M,'hopping_towards_fp.avi','fps',RATE/8,'compression','Cinepak');
end