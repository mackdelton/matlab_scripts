from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool


# initalize empty lists for data
a_waypoint_lat = []
a_waypoint_lon = []
b_waypoint_lat = []
b_waypoint_lon = []
a_lat = []
a_lon = []

# initialize dictionary to be used for hover data
data = {'lat': [], 'lon': []}
d1 = {i: [] for i in range(20)}
data.update(d1)

# read in data file and populate a list containing the contents of the file
data_file_in = open('pScheduleCalcBern_ucb_2017_07_25__13_19_09.txt.log', 'r')
lst = [line.strip() for line in data_file_in.readlines()]

# iterate over the list
for line in lst:
    # indication of A waypoints
    if 'Vehicle A Way Points' in line:
        for i in range(1, 21):
            # slice the line into lat and long based on the position of the = and , in the data. Append to lists.
            a_waypoint_lat.append(lst[lst.index(line) + i]
                                  [lst[lst.index(line) + i].index('=') + 2:lst[lst.index(line) + i].index(',')])
            a_waypoint_lon.append(lst[lst.index(line) + i]
                                  [lst[lst.index(line) + i].index(',') + 1:])
            
# commented out to use B waypoint average
        """
        elif 'Vehicle B Way Points' in line:
        for i in range(1, 21):
            b_waypoint_lat.append(lst[lst.index(line) + i]
                                  [lst[lst.index(line) + i].index('=') + 2:lst[lst.index(line) + i].index(',')])
            b_waypoint_lon.append(lst[lst.index(line) + i]
                                  [lst[lst.index(line) + i].index(',') + 1:])
        """
    
    # indication of B mean position
    elif 'Vehicle B Mean Position' in line:
        # slice the line into lat and lon based on the position of the = and , in the data. Append to lists.
        b_waypoint_lat.append(line[line.index('=') + 3:line.index(',')])
        b_waypoint_lon.append(line[line.index(',') + 1:])
    
    # indication of Vehicle A Position data
    elif 'Vehicle A Position' in line:
        # slice the line into lat and lon based on the position of the =, ,, and |. Append to lists.
        a_lat.append(line[line.index('=') + 1:line.index(',')])
        a_lon.append(line[line.index(',') + 1:line.index('|') - 1])
    """
    I'm working this part out right now. Basically it populates the dictionary of lists created earlier
    and adds the lines from the decision logic to plot during hover over. We need to discuss what you really
    want to see here as there are 274 of these events.
    # identify and format intelligence data for plotting
    elif line == 'UCBVec':
        for i in range(20):
            data[i].append(lst[lst.index(line) + i + 1])
        if not a_lat:
            data['lat'].append(0.0)
        else:
            data['lat'].append(a_lat[-1])
        if not a_lon:
            data['lon'].append(0.0)
        else:
            data['lon'].append(a_lat[-1])
    """
print(data)

# initialize html file for interactive bokeh plot
output_file("line.html")

p = figure(plot_width=800, plot_height=600)

# add a circle renderer with a size, color, and alpha
# Commented out to remove waypoints from plots
# p.circle(a_waypoint_lat, a_waypoint_lon, size=20, color="navy", alpha=0.5)
# p.circle(b_waypoint_lat, b_waypoint_lon, size=20, color="red", alpha=0.5)
p.circle(a_lat, a_lon, size=10, color='black', alpha=0.5)
# show the results
show(p)
