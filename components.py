import simulator, grn
import matplotlib.pyplot as plt

def add_rs_flip_flop(grn: grn.grn, prefix, r_name, s_name, kd=5, n=2, delta=0.2, alpha=10):
    q_name = f"{prefix}_Q"
    notQ_name = f"{prefix}_Q_"
    grn.add_species(q_name, delta)
    grn.add_species(notQ_name, delta)

    grn.add_gene(alpha, [
        {'name': r_name, 'type': -1, 'Kd': kd, 'n': n},
        {'name': notQ_name, 'type': -1, 'Kd': kd, 'n': n}
    ], [{'name': q_name}])

    grn.add_gene(alpha, [
        {'name': s_name, 'type': -1, 'Kd': kd, 'n': n},
        {'name': q_name, 'type': -1, 'Kd': kd, 'n': n}
    ], [{'name': notQ_name}])


def add_d_latch(grn:grn.grn, prefix, d_name, c_name, kd=5, n=2, delta=0.2, alpha=10):
    notD_name = f"{prefix}_notD"
    notC_name = f"{prefix}_notC"
    dandc_name = f"{prefix}_DandC"
    notdandc_name = f"{prefix}_notDandC"
    

    grn.add_species(notD_name, delta)
    grn.add_gene(alpha, [{'name': d_name, 'type': -1, 'Kd': kd, 'n': n}], [{'name': notD_name}])

    grn.add_species(notC_name, delta)
    grn.add_gene(alpha, [{'name': c_name, 'type': -1, 'Kd': kd, 'n': n}], [{'name': notC_name}])

    grn.add_species(dandc_name, delta)
    grn.add_gene(alpha, [
        {'name': notD_name, 'type': -1, 'Kd': kd, 'n': n},
        {'name': c_name, 'type': 1, 'Kd': kd, 'n': n}
    ], [{'name': dandc_name}])

    grn.add_species(notdandc_name, delta)
    grn.add_gene(alpha, [
        {'name': d_name, 'type': -1, 'Kd': kd, 'n': n},
        {'name': c_name, 'type': 1, 'Kd': kd, 'n': n}
    ], [{'name': notdandc_name}])

    add_rs_flip_flop(grn, prefix, notdandc_name, dandc_name, kd, n, delta, alpha)


def add_d_flip_flop(grn:grn.grn, prefix, d_name, c_name, kd=10, n=4, delta=0.7, alpha=40):
    dl1_name = f"{prefix}_DL1"
    q1_name = f"{dl1_name}_Q"

    add_d_latch(grn, dl1_name, d_name, c_name, kd, n, delta, alpha)

    d2_name = f"{prefix}_D2"
    c2_name = f"{prefix}_C2"
    
    grn.add_species(d2_name, delta)
    grn.add_species(c2_name, delta)

    grn.add_gene(alpha, [
        {'name': q1_name, 'type': 1, 'Kd': kd, 'n': n}
    ], [{'name': d2_name}])

    grn.add_gene(alpha, [
        {'name': f'{dl1_name}_notC', 'type': 1, 'Kd': kd, 'n': n}
    ], [{'name': c2_name}])

    dl2_name = f"{prefix}_DL2"
    add_d_latch(grn, dl2_name, d2_name, c2_name, kd, n, delta, alpha)
    



