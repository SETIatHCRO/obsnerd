#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# Author: wfarah
# GNU Radio version: 3.10.8.0

from PyQt5 import Qt
from gnuradio import qtgui
from datetime import datetime
from gnuradio import blocks
from gnuradio import fft
from gnuradio.fft import window
from gnuradio import gr
from gnuradio.filter import firdes
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import uhd
import time
import sip


def snipfcn_snippet_0(self):
    self.uhd_usrp_source.set_time_next_pps(uhd.time_spec(int(time.time())+1))
    self.uhd_usrp_source.set_lo_export_enabled(True, "lo1", 0)
    self.uhd_usrp_source.set_rx_lo_dist(True, "LO_OUT_0")
    self.uhd_usrp_source.set_rx_lo_dist(True, "LO_OUT_1")
    self.uhd_usrp_source.set_lo_source("external", "lo1", 0)
    self.uhd_usrp_source.set_lo_source("external", "lo1", 1)
    self.uhd_usrp_source.set_lo_source("external", "lo1", 2)
    self.uhd_usrp_source.set_lo_source("external", "lo1", 3)
    time.sleep(1)


def snippets_main_after_init(tb):
    snipfcn_snippet_0(tb)

class nrdz(gr.top_block, Qt.QWidget):

    def __init__(self, gaindB=30, samp_rate=122.88e6, src_name='unknown_source'):
        gr.top_block.__init__(self, "Not titled yet", catch_exceptions=True)
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Not titled yet")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except BaseException as exc:
            print(f"Qt GUI: Could not set Icon: {str(exc)}", file=sys.stderr)
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "nrdz")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Parameters
        ##################################################
        self.gaindB = gaindB
        self.samp_rate = samp_rate
        self.src_name = src_name

        ##################################################
        # Variables
        ##################################################
        self.tune = tune = 0
        self.nfft2 = nfft2 = 2**14
        self.nfft = nfft = 2**13
        self.fc = fc = 512e6

        ##################################################
        # Blocks
        ##################################################

        self.uhd_usrp_source = uhd.usrp_source(
            ",".join(("addr0=10.11.1.20, addr1=10.11.1.22, clock_source=external, time_source=external, ref_clk_freq=10e6", "")),
            uhd.stream_args(
                cpu_format="fc32",
                otw_format="sc16",
                args='',
                channels=list(range(0,4)),
            ),
        )
        self.uhd_usrp_source.set_samp_rate(samp_rate)
        self.uhd_usrp_source.set_time_unknown_pps(uhd.time_spec(0))

        self.uhd_usrp_source.set_center_freq(fc, 0)
        self.uhd_usrp_source.set_antenna('RX2', 0)
        self.uhd_usrp_source.set_gain(gaindB, 0)

        self.uhd_usrp_source.set_center_freq(fc, 1)
        self.uhd_usrp_source.set_antenna('RX2', 1)
        self.uhd_usrp_source.set_gain(gaindB, 1)

        self.uhd_usrp_source.set_center_freq(fc, 2)
        self.uhd_usrp_source.set_antenna('RX2', 2)
        self.uhd_usrp_source.set_gain(gaindB, 2)

        self.uhd_usrp_source.set_center_freq(fc, 3)
        self.uhd_usrp_source.set_antenna('RX2', 3)
        self.uhd_usrp_source.set_gain(gaindB, 3)
        self.qtgui_vector_sink_f_0 = qtgui.vector_sink_f(
            nfft,
            (-samp_rate/2+1.575e9),
            (samp_rate/nfft),
            "x-Axis",
            "y-Axis",
            "",
            1, # Number of inputs
            None # parent
        )
        self.qtgui_vector_sink_f_0.set_update_time(0.10)
        self.qtgui_vector_sink_f_0.set_y_axis((-20), 60)
        self.qtgui_vector_sink_f_0.enable_autoscale(False)
        self.qtgui_vector_sink_f_0.enable_grid(False)
        self.qtgui_vector_sink_f_0.set_x_axis_units("")
        self.qtgui_vector_sink_f_0.set_y_axis_units("")
        self.qtgui_vector_sink_f_0.set_ref_level(0)


        labels = ['', '', '', '', '',
            '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
            "magenta", "yellow", "dark red", "dark green", "dark blue"]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]

        for i in range(1):
            if len(labels[i]) == 0:
                self.qtgui_vector_sink_f_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_vector_sink_f_0.set_line_label(i, labels[i])
            self.qtgui_vector_sink_f_0.set_line_width(i, widths[i])
            self.qtgui_vector_sink_f_0.set_line_color(i, colors[i])
            self.qtgui_vector_sink_f_0.set_line_alpha(i, alphas[i])

        self._qtgui_vector_sink_f_0_win = sip.wrapinstance(self.qtgui_vector_sink_f_0.qwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_vector_sink_f_0_win)
        self.fft_vxx_0 = fft.fft_vcc(nfft, True, window.blackmanharris(nfft), True, 4)
        self.blocks_vector_to_streams_0 = blocks.vector_to_streams(gr.sizeof_float*nfft, 1)
        self.blocks_stream_to_vector_0 = blocks.stream_to_vector(gr.sizeof_gr_complex*1, nfft)
        self.blocks_null_sink_0 = blocks.null_sink(gr.sizeof_gr_complex*1)
        self.blocks_nlog10_ff_0 = blocks.nlog10_ff(10, nfft, 0)
        self.blocks_integrate_xx_0 = blocks.integrate_ff((int(0.2*samp_rate/nfft)), nfft)
        self.blocks_complex_to_mag_squared_0 = blocks.complex_to_mag_squared(nfft)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_complex_to_mag_squared_0, 0), (self.blocks_integrate_xx_0, 0))
        self.connect((self.blocks_integrate_xx_0, 0), (self.blocks_nlog10_ff_0, 0))
        self.connect((self.blocks_nlog10_ff_0, 0), (self.blocks_vector_to_streams_0, 0))
        self.connect((self.blocks_stream_to_vector_0, 0), (self.fft_vxx_0, 0))
        self.connect((self.blocks_vector_to_streams_0, 0), (self.qtgui_vector_sink_f_0, 0))
        self.connect((self.fft_vxx_0, 0), (self.blocks_complex_to_mag_squared_0, 0))
        self.connect((self.uhd_usrp_source, 3), (self.blocks_null_sink_0, 2))
        self.connect((self.uhd_usrp_source, 1), (self.blocks_null_sink_0, 1))
        self.connect((self.uhd_usrp_source, 2), (self.blocks_null_sink_0, 0))
        self.connect((self.uhd_usrp_source, 0), (self.blocks_stream_to_vector_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "nrdz")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_gaindB(self):
        return self.gaindB

    def set_gaindB(self, gaindB):
        self.gaindB = gaindB
        self.uhd_usrp_source.set_gain(self.gaindB, 0)
        self.uhd_usrp_source.set_gain(self.gaindB, 1)
        self.uhd_usrp_source.set_gain(self.gaindB, 2)
        self.uhd_usrp_source.set_gain(self.gaindB, 3)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.qtgui_vector_sink_f_0.set_x_axis((-self.samp_rate/2+1.575e9), (self.samp_rate/self.nfft))
        self.uhd_usrp_source.set_samp_rate(self.samp_rate)

    def get_src_name(self):
        return self.src_name

    def set_src_name(self, src_name):
        self.src_name = src_name

    def get_tune(self):
        return self.tune

    def set_tune(self, tune):
        self.tune = tune

    def get_nfft2(self):
        return self.nfft2

    def set_nfft2(self, nfft2):
        self.nfft2 = nfft2

    def get_nfft(self):
        return self.nfft

    def set_nfft(self, nfft):
        self.nfft = nfft
        self.qtgui_vector_sink_f_0.set_x_axis((-self.samp_rate/2+1.575e9), (self.samp_rate/self.nfft))

    def get_fc(self):
        return self.fc

    def set_fc(self, fc):
        self.fc = fc
        self.uhd_usrp_source.set_center_freq(self.fc, 0)
        self.uhd_usrp_source.set_center_freq(self.fc, 1)
        self.uhd_usrp_source.set_center_freq(self.fc, 2)
        self.uhd_usrp_source.set_center_freq(self.fc, 3)



def argument_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-g", "--gaindB", dest="gaindB", type=eng_float, default=eng_notation.num_to_str(float(30)),
        help="Set rx gain in decibels [default=%(default)r]")
    parser.add_argument(
        "-b", "--samp-rate", dest="samp_rate", type=eng_float, default=eng_notation.num_to_str(float(122.88e6)),
        help="Set sampling rate [Hz] [default=%(default)r]")
    parser.add_argument(
        "--src-name", dest="src_name", type=str, default='unknown_source',
        help="Set src_name [default=%(default)r]")
    return parser


def main(top_block_cls=nrdz, options=None):
    if options is None:
        options = argument_parser().parse_args()

    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls(gaindB=options.gaindB, samp_rate=options.samp_rate, src_name=options.src_name)
    snippets_main_after_init(tb)
    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    qapp.exec_()

if __name__ == '__main__':
    main()
