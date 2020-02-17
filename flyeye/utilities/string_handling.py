
def format_channel(channel):
    """ Returns string representation of <channel>. """

    if channel is None or channel == '':
        return None

    elif type(channel) == int:
        return 'ch{:d}'.format(channel)

    elif type(channel) != str:
        raise ValueError('Channel must be specified in string format.')

    elif 'ch' in channel:
        return channel

    elif channel.lower() in 'rgb':
        return format_channel('rgb'.index(channel.lower()))

    elif channel.lower() in ('red', 'green', 'blue'):
        return format_channel('rgb'.index(channel.lower()[0]))

    else:
        raise ValueError('Channel string not recognized.')


def standardize_channels(data):
    """
    Replace all channel names in <data> with standardized channel names.

    Args:

        data (pd.DataFrame)

    """
    for name in ('r', 'g', 'b', 'red', 'green', 'blue'):

        formatted_name = format_channel(name)

        if name in data.columns:
            data[formatted_name] = data[name]
            data.drop(name, axis=1, inplace=True)

        if name + '_std' in data.columns:
            data[formatted_name + '_std'] = data[name + '_std']
            data.drop(name + '_std', axis=1, inplace=True)

    return data
