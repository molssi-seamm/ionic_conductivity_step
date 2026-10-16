# -*- coding: utf-8 -*-

"""The graphical part of a IonicConductivity step"""

import pprint  # noqa: F401
import tkinter as tk
import tkinter.ttk as ttk

from .ionic_conductivity_parameters import IonicConductivityParameters
from seamm_md_analysis.step_support import run_own_md
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
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the IonicConductivity parameterss
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkIonicConductivity.reset_dialog
        """

        frame = super().create_dialog(title="Ionic Conductivity", widget="notebook")
        # make it large!
        screen_w = self.dialog.winfo_screenwidth()
        screen_h = self.dialog.winfo_screenheight()
        w = int(0.9 * screen_w)
        h = int(0.8 * screen_h)
        x = int(0.05 * screen_w / 2)
        y = int(0.1 * screen_h / 2)

        self.dialog.geometry(f"{w}x{h}+{x}+{y}")

        # Add a frame for the flowchart
        notebook = self["notebook"]
        flowchart_frame = ttk.Frame(notebook)
        self["flowchart frame"] = flowchart_frame
        notebook.add(flowchart_frame, text="Flowchart", sticky=tk.NSEW)

        self.tk_subflowchart = seamm.TkFlowchart(
            master=flowchart_frame,
            flowchart=self.node.subflowchart,
            namespace=self.namespace,
        )
        self.tk_subflowchart.draw()

        # Fill in the control parameters
        # Shortcut for parameters
        P = self.node.parameters

        # ionic_conductivity frame to isolate widgets
        frame = self["structure_frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Ionic Conductivity",
            labelanchor="n",
            padding=10,
        )

        for key in IonicConductivityParameters.parameters:
            self[key] = P[key].widget(frame)

        # and binding to change as needed
        for key in IonicConductivityParameters.parameters:
            if hasattr(self[key], "combobox"):
                self[key].combobox.bind(
                    "<<ComboboxSelected>>", self.reset_structure_frame
                )

        # and lay them out
        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        IonicConductivity parameters.

        This function simply lays them out row by row with
        aligned labels. You may wish a more complicated layout that
        is controlled by values of some of the control parameters.
        If so, edit or override this method

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkIonicConductivity.create_dialog
        """

        row = 0

        self["structure_frame"].grid(row=row, column=0, sticky=tk.EW, pady=10)
        row += 1
        self.reset_structure_frame()

        return row

    def reset_structure_frame(self, widget=None):
        """Layout the widgets in the ionic_conductivity frame
        as needed for the current state"""

        source = self["md source"].get()
        own_md = source == run_own_md

        # The sub-flowchart is used only when this step runs its own MD
        notebook = self["notebook"]
        if own_md:
            notebook.add(self["flowchart frame"])
        else:
            notebook.hide(self["flowchart frame"])

        frame = self["structure_frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        # Show the controls that apply, given the others
        P = self.node.parameters
        values = {
            key: self[key].get() for key in IonicConductivityParameters.parameters
        }
        row = 0
        widgets = []
        for key in IonicConductivityParameters.parameters:
            if P.applies(key, values):
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.W)
                widgets.append(self[key])
                row += 1

        sw.align_labels(widgets, sticky=tk.E)

        frame.columnconfigure(0, minsize=50)

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
