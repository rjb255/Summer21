# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# <h1>Drying Oils</h1>
# The fist exercise in Chemical Engineering: a simplified version of a drying process as given in fig 1.
# <br>
# <img src = "drying.png" width = 500px></img>
# <h2>Libraries</h2>

# %%
import numpy as np
import math
import random
import scipy
import csv
import pandas as pd
from pprint import pprint
from collections import deque

# %% [markdown]
# <h2>Classes<h2>

# %%


class UnitOp:
    def __init__(self: object, streams: list = [], nxt: list = [None], process=None, name: str = "", ticks: int = 1):
        self.feed = streams
        self.input = mixing(self.feed)
        self.next = nxt
        self.process = process
        self.__name__ = name
        self.output = [dict(self.input)]

    def measure(self, prop: str):
        if prop in self.input:
            if prop != "composition":
                error = np.random.normal() * self.input[prop]
                return self.input[prop] + error
            else:
                composition = {chem: comp * (1+np.random.normal())
                               for chem, comp in self.input[prop].items()}
                compTot = sum(composition.values())
                composition = {chem: comp/compTot for chem,
                               comp in composition.items()}
                return composition
        else:
            raise KeyError


class Reactor1(UnitOp):
    def __init__(self: object, streams: list = [], nxt: list = [None], name="", ticks: int = 0):
        super().__init__(streams=streams, nxt=nxt, process=self.reaction, name=name)

    def reaction(self):
        self.input = mixing(self.feed)
        output = dict(self.input)

        if self.input["composition"]["AA"] != 0 and self.input["composition"]["CO"] / self.input["composition"]["AA"] <= 2/3:
            change = self.input["composition"]["CO"] * 0.9
            output["composition"]["AA"] = self.input["composition"]["AA"] - change
            output["composition"]["ACO"] = self.input["composition"]["ACO"] + change
            output["composition"]["H2O"] = self.input["composition"]["H2O"] + change
            output["composition"]["CO"] = self.input["composition"]["CO"] - change
            self.output.append(output)
            self.next[0].feed = [self.output.pop(0)]
        else:
            print("Not enough AA")


class Reactor2(UnitOp):
    def __init__(self: object, streams: list = [], nxt: list = [None], size: int = 0, name: str = ""):
        super().__init__(streams=streams, nxt=nxt, process=self.react, name=name)
        if size == 0:
            self.data = pd.read_csv("smallReactor2Data.csv", delimiter=" ")
        else:
            self.data = pd.read_csv("largeReactor2Data.csv", delimiter=" ")

    def react(self):
        self.input = mixing(self.feed)
        self.output = dict(self.input)
        if 250+273.15 < self.input["temperature"] < 330+273.15:
            temperature = self.input["temperature"]
            conversion = np.interp(
                temperature, self.data["temperature"], self.data["conversion"])
            selectivity = np.interp(
                temperature, self.data["temperature"], self.data["selectivity"])

            molarIn = {chem: self.input["molarFlow"] * comp for chem,
                       comp in self.input["composition"].items()}
            molarOut = dict(molarIn)
            ACOreac = conversion*molarIn["ACO"]
            molarOut["ACO"] = molarIn["ACO"] - ACOreac
            molarOut["AA"] = molarIn["AA"] + ACOreac
            molarOut["DCO"] = molarIn["DCO"] + ACOreac
            dm = (molarOut["DCO"] - molarOut["Gum"]
                  * selectivity) / (selectivity + 2)
            molarOut["DCO"] = molarOut["DCO"] - 2*dm
            molarOut["Gum"] = molarOut["Gum"] + dm
            self.output["molarFlow"] = sum(molarOut.values())
            if self.output["molarFlow"] != 0:
                self.output["composition"] = {
                    chem: comp / self.output["molarFlow"] for chem, comp in molarOut.items()}
            self.output["molarMass"] = molarMass(self.output)

        else:
            print("warning: temp out of range, unknown behaviour")

        self.next[0].feed = [self.output]


class Filter(UnitOp):
    def __init__(self: object, streams: list = [], nxt: list = [None], name: str = ""):
        super().__init__(streams=streams, nxt=nxt, process=self.__filter, name=name)

    def __filter(self):
        self.input = mixing(self.feed)
        self.output = [dict(self.input), dict(self.input)]
        self.output[1]["molarFlow"] = self.input["composition"]["Gum"] * \
            self.input["molarFlow"]
        self.output[1]["composition"] = {
            c: 0 for c in self.output[1]["composition"].keys()}
        self.output[1]["composition"]["Gum"] = 1
        if self.output[1]["molarFlow"] != 0:
            self.output[1]["massFlow"] = self.input["MM"]["Gum"] * \
                self.output[1]["molarFlow"]
            self.output[1]["molarMass"] = self.output[1]["massFlow"] / \
                self.output[1]["molarFlow"]

        self.output[0]["molarFlow"] = self.input["molarFlow"] - \
            self.input["composition"]["Gum"] * self.input["molarFlow"]
        self.output[0]["massFlow"] = self.input["massFlow"] - \
            self.output[1]["massFlow"]
        self.output[0]["composition"]["Gum"] = 0
        if self.output[0]["molarFlow"] != 0:
            self.output[0]["molarMass"] = self.output[0]["massFlow"] / \
                self.output[0]["molarFlow"]
            self.output[0]["composition"] = {
                c: v/sum(self.output[0]["composition"].values()) for c, v in self.output[0]["composition"].items()}

        for i in range(2):
            self.next[i].feed = [self.output[i]]


class Column1(UnitOp):
    def __init__(self: object, streams: list = [], nxt: list = [None], name: str = ""):
        super().__init__(streams=streams, nxt=nxt, process=self.separate, name=name)

    def separate(self):
        self.input = mixing(self.feed)
        self.output = [dict(self.input), dict(self.input)]

        molarIn = {chem: self.input["molarFlow"] * comp for chem,
                   comp in self.input["composition"].items()}
        molar0 = dict(molarIn)
        molar0["AA"] = 0
        molar0["H2O"] = 0
        molar0["DCO"] *= 0.05
        molar0["CO"] *= 0.995
        molar0["ACO"] *= 0.998

        self.output[0]["molarFlow"] = sum(molar0.values())
        if sum(molar0.values()) != 0:
            self.output[0]["composition"] = {
                chem: mol / sum(molar0.values()) for chem, mol in molar0.items()}
            self.output[0]["molarMass"] = molarMass(self.output[0])
            self.output[0]["massFlow"] = self.output[0]["molarMass"] * \
                self.output[0]["molarFlow"]

        molar1 = {chem: molarIn[chem] - molar0[chem] for chem in molar0.keys()}
        self.output[1]["molarFlow"] = sum(molar1.values())
        if sum(molar1.values()) != 0:
            self.output[1]["composition"] = {
                chem: mol / sum(molar1.values()) for chem, mol in molar1.items()}
            self.output[1]["molarMass"] = molarMass(self.output[1])
            self.output[1]["massFlow"] = self.output[1]["molarMass"] * \
                self.output[1]["molarFlow"]

        self.next[0].feed[2] = self.output[0]
        self.next[1].feed[0] = self.output[1]


class Column2(UnitOp):
    def __init__(self: object, streams: list = [], nxt: list = [None], name: str = ""):
        super().__init__(streams=streams, nxt=nxt, process=self.separate, name=name)
        self.output = [dict(self.input), dict(self.input)]

    def separate(self):
        self.input = mixing(self.feed)
        self.output = [dict(self.input), dict(self.input)]

        molarIn = {chem: self.input["molarFlow"] * comp for chem,
                   comp in self.input["composition"].items()}
        molar0 = dict(molarIn)
        molar0["AA"] = 0
        molar0["H2O"] = 0
        self.output[0]["molarFlow"] = sum(molar0.values())
        if sum(molar0.values()) != 0:
            self.output[0]["composition"] = {
                chem: mol / sum(molar0.values()) for chem, mol in molar0.items()}
            self.output[0]["molarMass"] = molarMass(self.output[0])
            self.output[0]["massFlow"] = self.output[0]["molarMass"] / \
                self.output[0]["molarFlow"]

        molar1 = {chem: molarIn[chem] - molar0[chem] for chem in molar0.keys()}
        self.output[1]["molarFlow"] = sum(molar1.values())
        if sum(molar1.values()) != 0:
            self.output[1]["composition"] = {
                chem: mol / sum(molar1.values()) for chem, mol in molar1.items()}
            self.output[1]["molarMass"] = molarMass(self.output[1])
            self.output[1]["massFlow"] = self.output[1]["molarMass"] / \
                self.output[1]["molarFlow"]

        for i in range(2):
            self.next[i].feed = [self.output[i]]


class AcidSeparator(UnitOp):
    def __init__(self: object, streams: list = [], nxt: list = [None], name: str = ""):
        super().__init__(streams=streams, nxt=nxt, process=self.separate, name=name)

    def separate(self):
        self.input = mixing(self.feed)
        self.output = [dict(self.input), dict(self.input)]

        m1 = self.input["molarFlow"]
        x1 = self.input["composition"]["AA"]
        z1 = 0.08
        y1 = 0.99
        if y1 != z1:
            m2 = m1 * (x1 - z1) / (y1 - z1)
        m3 = m1 - m2

        self.output[0]["molarFlow"] = m2
        self.output[0]["composition"] = {
            c: 0 for c in self.output[0]["composition"].keys()}
        self.output[0]["composition"]["AA"] = y1
        self.output[0]["composition"]["H2O"] = 1 - y1
        self.output[0]["molarMass"] = molarMass(self.output[0])

        self.output[1]["molarFlow"] = m3
        self.output[1]["composition"] = {
            c: 0 for c in self.output[1]["composition"].keys()}
        self.output[1]["composition"]["AA"] = z1
        self.output[1]["composition"]["H2O"] = 1 - z1
        self.output[1]["molarMass"] = molarMass(self.output[1])

        self.next[0].feed[3] = self.output[0]
        self.next[1].feed[0] = self.output[1]


class Output(UnitOp):
    def __init__(self: object, streams: list = [], name: str = ""):
        super().__init__(streams=streams, name=name)


# %% [markdown]
# <h2>Functions</h2>

# %%
def mixing(streams):
    temperatures = np.array([s["temperature"] for s in streams])
    heatCapacity = np.array([s["Cp"] for s in streams])
    massFlow = np.array([s["massFlow"] for s in streams])
    molarFlowRate = np.array([s["molarFlow"] for s in streams])

    newProperties = dict(streams[0])

    if sum(massFlow) != 0:
        energy = sum(temperatures * heatCapacity * massFlow)
        compositions = {}
        for chem in streams[0]["composition"].keys():
            compositions[chem] = sum(
                [s["composition"][chem]*s["molarFlow"]/sum(molarFlowRate) for s in streams])
        newProperties["temperature"] = energy / sum(massFlow * heatCapacity)
        newProperties["Cp"] = sum(massFlow * heatCapacity) / sum(massFlow)
        newProperties["massFlow"] = sum(massFlow)
        newProperties["molarFlow"] = sum(molarFlowRate)
        newProperties["molarMass"] = sum(massFlow) / sum(molarFlowRate)
        newProperties["composition"] = compositions

    return newProperties


def molarMass(stream):
    return sum([comp*mm for comp, mm in zip(stream["MM"].values(), stream["composition"].values())])

# %% [markdown]
# <h2>Other</h2>


# %%
composition = {
    "CO": 0,
    "ACO": 0,
    "DCO": 0,
    "AA": 0,
    "H2O": 0,
    "Gum": 0,
}

MMs = {
    "CO":  312.5,
    "ACO": 354.5140,
    "DCO": 294.4620,
    "AA":  60.052,
    "H2O": 18.016,
    "Gum": 588.9240,
}

fluid = {
    "temperature": 0,
    "Cp": 0,
    "massFlow": 0,
    "molarFlow": 0,
    "molarMass": 0,
    "composition": None,
    "MM": MMs,
}

# %% [markdown]
# <h2>Creation</h2>

# %%
streams = [dict(fluid) for i in range(13)]
for s in streams:
    s["composition"] = dict(composition)

streams[0]["temperature"] = 573
streams[0]["Cp"] = 4200
streams[0]["molarFlow"] = 1660
streams[0]["composition"]["CO"] = 1
streams[0]["molarMass"] = molarMass(streams[0])
streams[0]["massFlow"] = streams[0]["molarFlow"] * streams[0]["molarMass"]

streams[1]["temperature"] = 573
streams[1]["Cp"] = 4200
streams[1]["molarFlow"] = 3660
streams[1]["composition"]["AA"] = 1
streams[1]["molarMass"] = molarMass(streams[1])
streams[1]["massFlow"] = streams[1]["molarFlow"] * streams[1]["molarMass"]


gum = Output([streams[5]], name="gum removal")
product = Output([streams[9]], name="product")
waste = Output([streams[12]], name="waste")
acid1 = AcidSeparator(streams=[streams[10]], nxt=[
                      None, waste], name="acid separator")
column2 = Column2(streams=[streams[8]], nxt=[product, acid1], name="column 2")
column1 = Column1(streams=[streams[6]], nxt=[None, column2], name="column 1")
filter1 = Filter(streams=[streams[4]], nxt=[gum, column1], name="filter")
reactor2 = Reactor2(streams=[streams[3]], nxt=[filter1], name="reactor 2")
reactor1 = Reactor1(streams=[streams[a] for a in [0, 1, 7, 11]], nxt=[
                    reactor2], name="reactor 1")
column1.next[0] = reactor1
acid1.next[0] = reactor1

objs = [acid1, column1, column2, filter1, reactor1, reactor2]
for i in range(4):
    for ob in objs:
        ob.process()
    print(column1.input["massFlow"] - column1.output[0]
          ["massFlow"] - column1.output[1]["massFlow"])


# %%
pprint(molarMass(reactor2.output))
pprint(reactor2.output)
