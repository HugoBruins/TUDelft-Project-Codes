"""
Python script to plot the data from https://github.com/HugoBruins/TUDelft-Project-Codes/blob/main/GMAT%20scripts/Project_Mars_Script_C07_Eclipselocator.script
It plots the sun exposure data from the first Eclipse onward. 
"""

import datetime as t
import matplotlib.pyplot as plt

# Changes the month as a word into a number
# Mar -> 03, Dec -> 12
def change_month(month: str) -> int:
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", 
              "Oct", "Nov", "Dec"]
    for n, index in enumerate(months):
        if month == index:
            return n+1
    return 0


# Changes the time format in the text file
# 07 Feb 2030 21:15:30.407 ---> 2030-2-7-21:15:30.407
def change_format(line: list) -> float:
    days = int(line[0])
    months = int(change_month(line[1]))
    years = int(line[2])
    time = line[3].split(':')
    hours = int(time[0])
    minutes = int(time[1])
    second = time[2].split('.')
    seconds = int(second[0])
    milliseconds = int(second[1]) * 1000
    total_time = t.datetime(years, months, days, hours, minutes, seconds, 
                            milliseconds)
    return total_time


# Read the file
filename = "C07EclipseData.txt"
file = open(filename, 'r')

# Format the file into a formatted begin_time, end_time and contact_time only
times = []
first = 0
for n,line in enumerate(file):
    # Not all lines contain data, so some have to be ignored
    splitline = line.split()
    if len(splitline) == 13:
        if first == 0:
            print("First eclipse: ", splitline[0:4])
            first = 1
        begin_time = splitline[0:4]
        end_time = splitline[4:8]
        begin_contact = change_format(begin_time)
        end_contact = change_format(end_time)
        contact_time = (end_contact - begin_contact).total_seconds()
        times.append([begin_contact, end_contact, contact_time])

max_length = len(times)
contact_time_array = [0]
total_time = [0]
for n,time in enumerate(times):
    if n != max_length-1:
        end_time_next = times[n+1][1]
        end_time_current = times[n][1]
        time_interval = (end_time_next - end_time_current).total_seconds()
        # Get the eclipse % and turn it into contact time and total time lists
        eclipse_time = time[2] / time_interval
        contact_time = 1 - eclipse_time
        contact_time *= 100
        # Make sure it shows up as sun synchronous, it has to do with how the
        # Data is set up in the file
        if 99.5 <= contact_time <= 100:
            contact_time_array[n] = 100
        contact_time_array.append(contact_time)
        # Add a total time array to plot on the x-axis
        total_time.append((total_time[n] + time_interval/(86400*365)))

del total_time[0]
del contact_time_array[0]
# Plotting the data
plt.plot(total_time, contact_time_array)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel("Earth Years", fontsize = 20)
plt.ylabel("Total sun exposure [%]",fontsize = 20)
plt.grid()
