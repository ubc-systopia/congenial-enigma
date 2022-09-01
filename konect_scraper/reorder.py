import config

def compute_random():
    return

def compute_rabbit():
    return

def compute_ordering(ord_str):
    match ord_str:
        case "rnd":
            compute_random()
        case "rbt":
            compute_rabbit()
        # case ""

def main(datasets, orde):
    settings = config.settings

    datasets_json_path = settings['datasets_json_path']
    graphs_dir = settings['graphs_dir']
    return

if __name__ == '__main__':
    main()
