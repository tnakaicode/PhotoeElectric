import numpy as np
import scipy.constants as cnt

SI = {
    'nm': 10**(-9), 'um': 10**(-6), 'mm': 0.001, 'cm': 0.01, 'm': 1.0, 'km': 1000.,
    'THz': 10**12, 'GHz': 10**9, 'MHz': 10**6, 'kHz': 10**3, 'Hz': 1.0,
    'ns': 10**(-9), 'us': 10**(-6), 'ms': 10**(-3), 's': 1.0,
    'min': 60.0, 'h': 60.0 * 60.0, 'day': 60.0 * 60.0 * 24.0,
    'month': 60.0 * 60.0 * 24.0 * 30,
    'year': 60.0 * 60.0 * 24.0 * 365
}


def convert_SI(val_unit, unit_in="cm", unit_out="mm"):
    val = float(val_unit.split(unit_in)[0])
    return val * SI[unit_in] / SI[unit_out]


def convert(unit_in="cm", unit_out="mm"):
    return SI[unit_in] / SI[unit_out]


def convert_freq_to_wave(num=100.0, unit_in="GHz", unit_out="mm"):
    """
    v = f * lambda
    """
    freq = num * SI[unit_in] / SI["Hz"]
    wave = cnt.c / freq
    return wave * SI["m"] / SI[unit_out]


def convert_wave_to_freq(num=1.0, unit_in="mm", unit_out="GHz"):
    """
    v = f * lambda
    """
    wave = num * SI[unit_in] / SI["m"]
    freq = cnt.c / wave
    return freq * SI["Hz"] / SI[unit_out]


def convert_freq_to_time(num=1.0, unit_in="Hz", unit_out="s"):
    freq = num * SI[unit_in] / SI["Hz"]
    time_s = 1 / freq
    return time_s * SI["s"] / SI[unit_out]


def convert_time_to_freq(num=1.0, unit_in="s", unit_out="Hz"):
    time_s = num * SI[unit_in] / SI["s"]
    freq = 1 / time_s
    return freq * SI["Hz"] / SI[unit_out]


if __name__ == "__main__":
    # freq -> time
    print(convert_freq_to_time(1.0, "Hz", "s"), "s")
    print(convert_freq_to_time(2.0, "Hz", "s"), "s")
    print(convert_freq_to_time(2.0, "kHz", "ms"), "ms")
    print(convert_freq_to_time(10.0, "kHz", "us"), "us")
    # time -> freq
    print(convert_time_to_freq(1.0, "s", "Hz"), "Hz")
    print(convert_time_to_freq(0.5, "s", "Hz"), "Hz")
    print(convert_time_to_freq(100, "ms", "Hz"), "Hz")
    print(convert_time_to_freq(100, "us", "kHz"), "kHz")
    print(convert_time_to_freq(100, "ns", "MHz"), "MHz")
