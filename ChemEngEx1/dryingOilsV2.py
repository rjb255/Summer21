# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
from IPython import get_ipython

# %% [markdown]
# <h1>Drying Oils</h1>
# The first exercise in Chemical Engineering: a simplified version of a drying process as given in fig 1.
# <br>
# <img src = "drying.png" width = 500px></img>
# <h2>Libraries</h2>

# %%
from abc import abstractmethod, ABC
from pprint import pprint
from collections import deque

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn import metrics
from scipy import optimize as opt

# %% [markdown]
# <h2>Units<h2>

# %%

class UnitOp(ABC):
    # Generic Unit Operation
    # Each unit has an input stream, outputs, following units, a name, and a process function
    # A unit encompasses the process unit and the stream leaving up until the following unit operations
    # The input of each stream can be measured with an uncertainty attached
    # The process each unit carries out is given in the attached documents

    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
    ):
        self.inputObjects = [None] * len(inputShapes)
        self.outputObjects = [None] * len(outputShapes)
        self.__meta__ = meta
        self.input = [[None]*a for a in inputShapes]
        self.output = [[None]*a for a in outputShapes]
        self.profitModel = lambda x: 0
        self.modelProcess = lambda x: 0

    def addInput(self, objFrom, addition):
        for index, obj in enumerate(self.inputObjects):
            if objFrom is obj:
                break
        self.input.insert(addition)

    @abstractmethod
    def process(self):
        for _input, _output in zip(self.input, self.output):
            _output.insert(0, _input.pop())
        for outObj, _output in zip(self.outputObjects, self.output):
            outObj.addInput(self, _output.pop())

    
    @abstractmethod
    def profit(self):
        pass

    # Needs Modification
    def measure(self, prop: str):
        # Return the value of the specified property with an error attached
        if prop in self.input:
            if prop != "composition" and prop != "MM":
                error = np.random.normal() * self.input[prop]
                return self.input[prop] + error/50
            elif prop == "MM":
                return self.input[prop]
            else:
                composition = {
                    chem: comp * (1 + np.random.normal())
                    for chem, comp in self.input[prop].items()
                }
                compTot = sum(composition.values())
                if compTot == 0:
                    compTot = float("inf")
                composition = {
                    chem: comp / compTot for chem, comp in composition.items()
                }
                return composition
        else:
            raise KeyError


class Pipe(UnitOp):
    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
    ):
        super().__init__(inputShapes, outputShapes, meta)

    def profit(self):
        return 0

    def process(self):
        super().process()


class Reactor1(UnitOp):
    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
    ):
        super().__init__(inputShapes, outputShapes, meta)


    def process(self):
        # Parameters
        X = 0.9 # Conversion
        feed = mixing([inp.pop() for inp in self.input])
        output = dict(feed)
        if feed is None:
            self.output[0].addInput(self, output)
        elif (
            feed["composition"]["CO"] * 1.5 < feed["composition"]["AA"]
        ):
            # Reaction1: AA + CO -> ACO + H2O with 90% conversion of CO and 50% excess AA
            change = self.input["composition"]["CO"] * X
            output["composition"]["AA"] = self.input["composition"]["AA"] - change
            output["composition"]["ACO"] = self.input["composition"]["ACO"] + change
            output["composition"]["H2O"] = self.input["composition"]["H2O"] + change
            output["composition"]["CO"] = self.input["composition"]["CO"] - change
        else:
            raise Warning("Not enough AA")
        self.output[0].addInput(self, output)

    def profit(self):
        return (
            -self.input["massFlow"] *
            (self.input["temperature"] - 293.15) / 10000
            - self.feed[2]["massFlow"] / 1000
            - self.feed[3]["massFlow"] / 1000
        )


class Reactor2(UnitOp):
    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
        size = 0
    ):
        super().__init__(inputShapes, outputShapes, meta)
        if size == 0:
            self.data = pd.read_csv("smallReactor2Data.csv", delimiter=" ")
        else:
            self.data = pd.read_csv("largeReactor2Data.csv", delimiter=" ")

    #Change Interpolation and Range stuff
    def process(self):
        # Parameters
        lowerWT = 270 + 273.15
        upperWT = 310 + 273.15
        feed = mixing([inp.pop() for inp in self.input])
        output = dict(feed)
        temperature = output["temperature"]
        if 250 + 273.15 < temperature < 330 + 273.15:
            if lowerWT > temperature:
                raise Warning("TEMPERATURE LOW")
            elif upperWT < temperature:
                raise Warning("TEMPERATURE HIGH")

            # From empirical data
            conversion = np.interp(
                temperature, self.data["temperature"], self.data["conversion"]
            )
            selectivity = np.interp(
                temperature, self.data["temperature"], self.data["selectivity"]
            )
            # Molar flow rates of each species
            molarIn = {
                chem: output["molarFlow"] * comp
                for chem, comp in output["composition"].items()
            }
            molarOut = dict(molarIn)

            # Reaction1: ACO -> DCO + AA
            ACOreac = conversion * molarIn["ACO"]
            molarOut["ACO"] -= ACOreac
            molarOut["AA"] += ACOreac
            molarOut["DCO"] += ACOreac

            # Reaction2: 2DCO -> Gum
            dm = (molarOut["DCO"] - molarOut["Gum"]
                  * selectivity) / (selectivity + 2)
            molarOut["DCO"] = molarOut["DCO"] - 2 * dm
            molarOut["Gum"] = molarOut["Gum"] + dm
            output["molarFlow"] = sum(molarOut.values())

            # Normalisation of molar flows
            if output["molarFlow"] != 0:
                output["composition"] = {
                    chem: comp / output["molarFlow"] for chem, comp in molarOut.items()
                }
            output["molarMass"] = molarMass(output)
        else:
            raise Exception("TEMPERATURE")
        self.output[0].addInput(self, output)

    def profit(self):
        return -self.input["massFlow"] * (self.input["temperature"] - 293.15) / 10000


class Filter(UnitOp):
    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
    ):
        super().__init__(inputShapes, outputShapes, meta)

    def process(self):
        self.input = mixing(self.feed)
        output = [dict(self.input), dict(self.input)]

        output[0]["molarFlow"] *= self.input["composition"]["Gum"]
        output[0]["composition"] = {
            c: 0 for c in output[0]["composition"].keys()}
        output[0]["composition"]["Gum"] = 1
        output[0]["massFlow"] = self.input["MM"]["Gum"] * \
            output[0]["molarFlow"]
        output[0]["molarMass"] = molarMass(output[0])

        output[1]["molarFlow"] = self.input["molarFlow"] - \
            output[0]["molarFlow"]
        output[1]["massFlow"] = self.input["massFlow"] - output[0]["massFlow"]
        output[1]["composition"]["Gum"] = 0
        output[1]["molarMass"] = molarMass(output[1])
        if output[1]["molarFlow"] != 0:
            output[1]["composition"] = {
                c: v / sum(output[1]["composition"].values())
                for c, v in output[1]["composition"].items()
            }
        self.output.append(output)
        for i, out in enumerate(self.output.pop(0)):
            self.next[i].feed = [out]

    def profit(self):
        return -self.output[-1][0]["massFlow"] * 40


class Column1(UnitOp):
    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
    ):
        super().__init__(inputShapes, outputShapes, meta)

    def separate(self):
        self.input = mixing(self.feed)
        output = [dict(self.input), dict(self.input)]

        molarIn = {
            chem: output[0]["molarFlow"] * comp
            for chem, comp in output[0]["composition"].items()
        }
        # Heavy Components (0.05 DCO, 0.995 CO, 0.998 ACO)
        molar0 = dict(molarIn)
        molar0["AA"] = 0
        molar0["H2O"] = 0
        molar0["DCO"] *= 0.05
        molar0["CO"] *= 0.995
        molar0["ACO"] *= 0.998

        output[0]["molarFlow"] = sum(molar0.values())
        if output[0]["molarFlow"] != 0:
            output[0]["composition"] = {
                chem: mol / output[0]["molarFlow"] for chem, mol in molar0.items()
            }
        output[0]["molarMass"] = molarMass(output[0])
        output[0]["massFlow"] = output[0]["molarMass"] * output[0]["molarFlow"]

        # Top Components
        molar1 = {chem: molarIn[chem] - molar0[chem] for chem in molar0.keys()}
        output[1]["molarFlow"] = sum(molar1.values())
        if output[1]["molarFlow"] != 0:
            output[1]["composition"] = {
                chem: mol / output[1]["molarFlow"] for chem, mol in molar1.items()
            }
        output[1]["molarMass"] = molarMass(output[1])
        output[1]["massFlow"] = output[1]["molarMass"] * output[1]["molarFlow"]

        self.output.append(output)
        out = self.output.pop(0)
        self.next[0].feed[2] = out[0]
        self.next[1].feed[0] = out[1]

    def profit(self):
        return -self.input["massFlow"] / 25


class Column2(UnitOp):
    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
    ):
        super().__init__(inputShapes, outputShapes, meta)

    def separate(self):
        self.input = mixing(self.feed)
        output = [dict(self.input), dict(self.input)]

        molarIn = {
            chem: output[0]["molarFlow"] * comp
            for chem, comp in output[0]["composition"].items()
        }
        # Bottom fraction: No AA or H2O
        molar0 = dict(molarIn)
        molar0["AA"] = 0
        molar0["H2O"] = 0
        output[0]["molarFlow"] = sum(molar0.values())
        if output[0]["molarFlow"] != 0:
            output[0]["composition"] = {
                chem: mol / output[0]["molarFlow"] for chem, mol in molar0.items()
            }
        output[0]["molarMass"] = molarMass(output[0])
        output[0]["massFlow"] = output[0]["molarMass"] * output[0]["molarFlow"]

        # Top Fraction
        molar1 = {chem: molarIn[chem] - molar0[chem] for chem in molar0.keys()}
        output[1]["molarFlow"] -= output[0]["molarFlow"]
        if output[1]["molarFlow"] != 0:
            output[1]["composition"] = {
                chem: mol / output[1]["molarFlow"] for chem, mol in molar1.items()
            }
        output[1]["molarMass"] = molarMass(output[1])
        output[1]["massFlow"] -= output[0]["massFlow"]
        self.output.append(output)

        for i, out in enumerate(self.output.pop(0)):
            self.next[i].feed = [out]

    def profit(self):
        return -self.input["massFlow"] / 100


class AcidSeparator(UnitOp):
    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
    ):
        super().__init__(inputShapes, outputShapes, meta)

    def separate(self):
        self.input = mixing(self.feed)
        output = [dict(self.input), dict(self.input)]

        m1 = self.input["molarFlow"]
        x1 = self.input["composition"]["AA"]
        z1 = 0.08
        y1 = 0.99
        if y1 != z1:
            m2 = m1 * (x1 - z1) / (y1 - z1)
        m3 = m1 - m2

        output[0]["molarFlow"] = m2
        output[0]["composition"] = {
            c: 0 for c in output[0]["composition"].keys()}
        output[0]["composition"]["AA"] = y1
        output[0]["composition"]["H2O"] = 1 - y1
        output[0]["molarMass"] = molarMass(output[0])

        output[1]["molarFlow"] = m3
        output[1]["composition"] = {
            c: 0 for c in output[1]["composition"].keys()}
        output[1]["composition"]["AA"] = z1
        output[1]["composition"]["H2O"] = 1 - z1
        output[1]["molarMass"] = molarMass(output[1])

        self.output.append(output)
        out = self.output.pop(0)

        self.next[0].feed[3] = out[0]
        self.next[1].feed[0] = out[1]

    def profit(self):
        return -self.input["massFlow"] * 0.15


class Output(UnitOp):
    def __init__(
        self: object,
        /,
        inputShapes: list = [1],
        outputShapes: list = [1],
        meta = None,
    ):
        super().__init__(inputShapes, outputShapes, meta)
    
    def process():
        pass

    def profit(self):
        if self.__name__ == "product":
            return (
                self.feed[0]["composition"]["DCO"]
                * self.feed[0]["molarFlow"]
                * self.feed[0]["MM"]["DCO"]
            )
        else:
            return 0

# %% [markdown]
# <h2>Functions</h2>

# %%


def mixing(streams):
    temperatures = np.array([s["temperature"] for s in streams])
    heatCapacity = np.array([s["Cp"] for s in streams])
    massFlow = np.array([s["massFlow"] for s in streams])
    mm = np.array([molarMass(s) for s in streams])
    molarFlowRate = massFlow/mm
    molarFlowRate[np.isnan(molarFlowRate)] = 0
    for s, mf in zip(streams, molarFlowRate):
        s["molarFlow"] = mf
    newProperties = dict(streams[0])

    if sum(massFlow) != 0:
        energy = sum(temperatures * heatCapacity * massFlow)
        compositions = {}
        for chem in streams[0]["composition"].keys():
            compositions[chem] = sum(
                [
                    s["composition"][chem] *
                    s["molarFlow"] / sum(molarFlowRate)
                    for s in streams
                ]
            )
        newProperties["temperature"] = energy / sum(massFlow * heatCapacity)
        newProperties["Cp"] = sum(massFlow * heatCapacity) / sum(massFlow)
        newProperties["massFlow"] = sum(massFlow)
        newProperties["molarFlow"] = sum(molarFlowRate)
        newProperties["molarMass"] = sum(massFlow) / sum(molarFlowRate)
        newProperties["composition"] = compositions

    return newProperties


def molarMass(stream):
    return sum(
        [
            comp * mm
            for mm, comp in zip(stream["MM"].values(), stream["composition"].values())
        ]
    )


def grossProfit(streams, units):
    profit = 0
    for s in streams:
        profit -= s["massFlow"] * s["composition"]["AA"] * 0.9
        profit += s["massFlow"] * s["composition"]["CO"] * 1.2
    # output = lambda u: u.input["composition"]["DCO"] * 2
    for unit in units:
        profit += unit.profit()
    return profit

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
    "CO": 312.5,
    "ACO": 354.5140,
    "DCO": 294.4620,
    "AA": 60.052,
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


class Model1:

    def __init__(self):
        self.streams = [dict(fluid) for i in range(13)]
        for s in self.streams:
            s["composition"] = dict(composition)

        self.streams[0]["temperature"] = 573
        self.streams[0]["Cp"] = 4200
        self.streams[0]["molarFlow"] = 1660
        self.streams[0]["composition"]["CO"] = 1
        self.streams[0]["molarMass"] = molarMass(self.streams[0])
        self.streams[0]["massFlow"] = self.streams[0]["molarFlow"] * \
            self.streams[0]["molarMass"]

        self.streams[1]["temperature"] = 573
        self.streams[1]["Cp"] = 4200
        self.streams[1]["molarFlow"] = 3660
        self.streams[1]["composition"]["AA"] = 1
        self.streams[1]["molarMass"] = molarMass(self.streams[1])
        self.streams[1]["massFlow"] = self.streams[1]["molarFlow"] * \
            self.streams[1]["molarMass"]

        self.gum = Output([self.streams[5]], name="gum removal")
        self.product = Output([self.streams[9]], name="product")
        self.waste = Output([self.streams[12]], name="waste")
        self.acid1 = AcidSeparator(streams=[self.streams[10]], nxt=[
            None, self.waste], name="acid separator")
        self.column2 = Column2(streams=[self.streams[8]], nxt=[
                               self.product, self.acid1], name="column 2")
        self.column1 = Column1(streams=[self.streams[6]], nxt=[
                               None, self.column2], name="column 1", ticks=3)
        self.filter1 = Filter(streams=[self.streams[4]], nxt=[
                              self.gum, self.column1], name="filter")
        self.reactor2 = Reactor2(streams=[self.streams[3]], nxt=[
                                 self.filter1], name="reactor 2")
        self.reactor1 = Reactor1(
            streams=[self.streams[a] for a in [0, 1, 7, 11]], nxt=[self.reactor2], name="reactor 1"
        )
        self.column1.next[0] = self.reactor1
        self.acid1.next[0] = self.reactor1

    def iterate(self):
        objs = [
            self.acid1,
            self.column1,
            self.filter1,
            self.reactor1,
            self.reactor2,
        ]
        for ob in objs:
            ob.process()


ex1 = Model1()

# %% [markdown]
# <h2> Controller </h2>

# %%


class Controller:
    def __init__(self, *, model, goal):
        self.objs = lambda model: [
            model.gum,
            model.product,
            model.waste,
            model.acid1,
            model.column1,
            model.filter1,
            model.reactor1,
            model.reactor2,
        ]
        self.model = model
        self.adjust = self.actions(model)
        self.profit = goal
        self.goal = self.target(model)
        self.prop = self.props(model)
        self.his = []

    def actions(self, model):
        adjustments = [[ex1.streams[1], ["massFlow", "temperature"]]]
        return adjustments

    def props(self, model):
        props = [
            [model.streams[0], [m for m in model.streams[0].keys() if m !=
                                "MM" and m != "composition"]],
            [model.streams[1], [m for m in model.streams[1].keys() if m !=
                                "MM" and m != "composition"]],
        ]
        props += [[obj, [m for m in obj.input.keys() if m != "MM" and m != "composition"]]
                  for obj in self.objs(model)]
        return props

    def target(self, model):
        return self.profit(model.streams[0:2], self.objs(model))

    def array(self, arr, index):
        if isinstance(index, int):
            for p in arr:
                index -= len(p[1])
                if index < 0:
                    break
            if isinstance(p[0], dict):
                return p[0][p[1][index + len(p[1])]]
            else:
                return p[0].input[p[1][index + len(p[1])]]
        else:
            raise IndexError

    def length(self, arr):
        length = 0
        for p in arr:
            length += len(p[1])
        return length

    def start(self, iterations):
        for i in range(iterations):
            control.adjust[0][0][control.adjust[0][1][0]
                                 ] = control.prop[8][0].measure(control.prop[8][1][2])
            self.model.iterate()
            measurements = [[], [], []]
            for pro in self.prop:
                for aspect in pro[1]:
                    if isinstance(pro[0], dict):
                        measurements[0].append(pro[0][aspect])
                    else:
                        measurements[0].append(pro[0].measure(aspect))
            measurements[1].append(control.adjust[0][0]
                                   [control.adjust[0][1][0]])
            measurements[2].append(self.target(self.model))
            print(round(measurements[2][-1], -3))
            self.his.append(measurements)
        print("\n")
        self.con(2)

    def con(self, iteration):
        for i in range(iteration):
            self.model.iterate()
            X = [h[0] + h[1] for h in self.his]
            Y = [h[2][0] for h in self.his]

            model = RandomForestRegressor(n_estimators=10)
            pip = Pipeline(steps=[('model', model)])
            currentModel = pip.fit(X, Y)
            currentMeasurements = []
            for pro in self.prop:
                for aspect in pro[1]:
                    if isinstance(pro[0], dict):
                        currentMeasurements.append(pro[0][aspect])
                    else:
                        currentMeasurements.append(pro[0].measure(aspect))
            adjust = control.adjust[0][0][control.adjust[0][1][0]]

            def toMin(x, xi, model):
                X = x + xi
                X = np.reshape(X, (1, -1))
                return -model.predict(X)
            adjust = opt.fmin(lambda x: toMin(
                list(x), currentMeasurements, currentModel), adjust, disp=0)
            control.adjust[0][0][control.adjust[0][1][0]] = adjust
            print(self.target(self.model))


control = Controller(goal=profit, model=ex1)


# %%
control.start(20)


# %%
Q = np.zeros(shape=[4, 2])
control.length(control.prop)
