# -*- coding: utf-8 -*-

"""The graphical part of a IonicConductivity step"""

import pprint  # noqa: F401
import tkinter as tk

from .ionic_conductivity_parameters import IonicConductivityParameters
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkIonicConductivity(seamm.TkNode):
    """
    The graphical part of a IonicConductivity step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in IonicConductivity parameterss.py

    See Also
    --------
    IonicConductivity, TkIonicConductivity,
    IonicConductivityParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        namespace="org.molssi.seamm.tk",
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.namespace = namespace
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )
        self.create_dialog()

    def create_dialog(self):
        """Create the dialog, with the widgets for the parameters.

        See Also
        --------
        TkIonicConductivity.reset_dialog
        """
        frame = super().create_dialog(title="Ionic Conductivity")

        P = self.node.parameters
        for key in IonicConductivityParameters.parameters:
            self[key] = P[key].widget(frame)
            if hasattr(self[key], "combobox"):
                self[key].combobox.bind("<<ComboboxSelected>>", self.reset_dialog)

        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Lay out the widgets that apply, given the choices made.

        Parameters
        ----------
        widget : Tk Widget = None
        """
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        P = self.node.parameters
        values = {
            key: self[key].get() for key in IonicConductivityParameters.parameters
        }
        row = 0
        widgets = []
        for key in IonicConductivityParameters.parameters:
            if P.applies(key, values):
                self[key].grid(row=row, column=0, sticky=tk.W)
                widgets.append(self[key])
                row += 1
        sw.align_labels(widgets, sticky=tk.E)
        return row

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkIonicConductivity.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit...", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
