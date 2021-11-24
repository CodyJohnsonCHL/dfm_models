"""Code for working with OpenDA

cody.l.johnson@erdc.dren.mil

"""


def writeNoosTs(data, Location, Unit, fn):
    """write time series data in OpenDA ExchangeObject format (NOOS)

    data = pd.Series with datatime axis and value

    """

    with open(fn, "w") as f:
        f.write("#------------------------------------------------------\n")
        f.write("#------------------------------------------------------\n")
        f.write(f"# Location : {Location}\n")
        f.write("# Position : (0.0,0.0)\n")
        f.write("# Source : null\n")
        f.write(f"# Unit : {Unit}\n")
        f.write("# Analyse time : null\n")
        f.write("# Timezone : null\n")
        f.write("#------------------------------------------------------\n")

        for time, value in data.iteritems():
            strftime = time.strftime("%Y%m%d%H%M")
            f.write(f"{strftime}\t{value}\n")
