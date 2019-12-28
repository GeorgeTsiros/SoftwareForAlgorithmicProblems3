import csv
import sys, getopt

predictions = []
inputfile = sys.argv[2]
with open(inputfile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:        
        line_count += 1
        predictions.append(row)
    print(predictions)
    print(f'Processed {line_count} lines.')
    print(len(predictions))