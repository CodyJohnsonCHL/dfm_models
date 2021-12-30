"""Code for working with OpenDA

cody.l.johnson@erdc.dren.mil

"""

from datetime import datetime as dt


def writeNoosTs(data, Location, x, y, Unit, fn):
    """write time series data in OpenDA ExchangeObject format (NOOS)

    data = pd.Series with datatime axis and value

    """
    now = dt.now().strftime("%Y-%m-%d %H:%M")

    with open(fn, "w") as f:
        f.write("#======================================================\n")
        f.write(f"# Generated on {now} \n")
        f.write("#======================================================\n")
        f.write(f"# Location         : {Location}\n")
        f.write(f"# Position         : ({x:.05f},{y:.05f})\n")
        f.write("# Source           : observed\n")
        f.write(f"# Unit             : {Unit}\n")
        f.write("# Analyse time     : null\n")
        f.write("#======================================================\n")

        for time, value in data.iteritems():
            strftime = time.strftime("%Y%m%d%H%M")
            f.write(f"{strftime}\t{value}\n")


def createNoosConfigFile(stdevs, stationNames, fn):

    with open(fn, "w") as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write(
            '<noosObserver xmlns="http://www.openda.org" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.openda.org http://schemas.openda.org/observation/noosObservations.xsd">\n'
        )
        f.write("\n")

        for stdev, stationName in zip(stdevs, stationNames):
            writeNoosObserver(f, stdev, stationName)

        f.write("\n")
        f.write("</noosObserver>")


def writeNoosObserver(f, stdev, noosFn):
    f.write(f'\t<timeSeries status="use" standardDeviation="{stdev:.02f}">\n')
    f.write(f"\t\t{noosFn}\n")
    f.write(f"\t</timeSeries>\n")
