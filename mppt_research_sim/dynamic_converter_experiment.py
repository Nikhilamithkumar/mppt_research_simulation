"""Optional averaged-boost-converter check using the same MPPT references."""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pv_model import PVModule, make_string_iv_curve
from mppt_algorithms import make_algorithm
from boost_converter import AveragedBoostConverter, BoostParameters

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"; OUT.mkdir(exist_ok=True)
MODULE=PVModule()


def run_dynamic(name, alg_name, n_samples=120, change_at=60, seed=2, partial=False):
    # This function runs the same MPPT logic against an averaged boost-converter
    # plant rather than the ideal instantaneous voltage-actuation benchmark. The
    # dynamic states make the plant behave more like a real converter and reveal
    # controller/plant interactions that the static benchmark would hide.
    params=BoostParameters()
    profiles={}
    profile_mpps={}
    def prof(k):
        if partial: return (1000.0,1000.0,300.0)
        g=1000.0 if k<change_at else 600.0
        return (g,g,g)
    def curve(gs):
        if gs not in profiles:
            profiles[gs]=make_string_iv_curve(gs,25.0,MODULE,bypass_v=0.7,points=4000)
            profile_mpps[gs]=profiles[gs].mpp()["power_w"]
        return profiles[gs]
    c0=curve(prof(0))
    vmin=0.0; vmax=0.97*c0.voc_v
    if partial:
        start_v=0.92*c0.voc_v
    else:
        start_v=0.55*c0.voc_v
    alg=make_algorithm(alg_name,vmin,vmax,start_v,po_step=0.7,
                       inc_deadband=0.005,pso_particles=6,pso_iterations=8,seed=seed)
    vref=alg.initial_voltage()
    converter=AveragedBoostConverter(c0,vref,params=params,dt_s=100e-6)
    dt_mppt_s=0.05
    power=[]; mpp=[]; vpv_hist=[]; vout_hist=[]; duty_hist=[]
    for k in range(n_samples):
        # One MPPT update cycle: sample the converter states, compute the current
        # power, ask the MPPT routine for the next voltage reference, and then let
        # the averaged converter settle under that command for a brief interval.
        gs=prof(k); c=curve(gs)
        vpv,ipv,pin,vout=converter.measured(c)
        power.append(pin); mpp.append(profile_mpps[gs]); vpv_hist.append(vpv);vout_hist.append(vout)
        vref=float(alg.next_voltage(vpv,ipv,pin))
        converter.advance(vref,c,dt_mppt_s)
        duty_hist.append(converter.last_duty or 0.0)
    power=np.array(power); mpp=np.array(mpp); vpv_hist=np.array(vpv_hist);vout_hist=np.array(vout_hist);duty_hist=np.array(duty_hist)
    return {"name":name,"algorithm":alg_name,"power":power,"mpp":mpp,"vpv":vpv_hist,
            "vout":vout_hist,"duty":duty_hist,"efficiency":float(power.sum()/mpp.sum()*100),
            "sample_period_s":dt_mppt_s}


def main():
    # Run the optional averaged-converter check for both an irradiance step and a
    # partial-shading case. This is a separate validation mode from the ideal
    # benchmark and is useful when checking how controller decisions interact with
    # converter dynamics.
    start=time.time(); runs=[]
    for condition,partial in [("irradiance_step",False),("partial_shading",True)]:
        for alg in ["P&O","Incremental Conductance","PSO"]:
            run=run_dynamic(condition,alg,partial=partial,seed=5)
            runs.append(run)
            print(condition,alg,"efficiency",round(run['efficiency'],2),"final power",round(float(run['power'][-1]),2),"final MPP",round(float(run['mpp'][-1]),2))
    summary={"model":"averaged CCM boost converter; ideal switch/diode with small inductor resistance",
             "outer_mppt_sample_period_s":0.05,"results":[],"runtime_seconds":None,
             "note":"This optional dynamic test is a numerical model, not hardware validation. Verify controller settings before using results as research evidence."}
    fig,axs=plt.subplots(2,2,figsize=(12,8),sharex="col")
    for r in runs:
        idx=0 if r['name']=="irradiance_step" else 1
        color={"P&O":"tab:blue","Incremental Conductance":"tab:orange","PSO":"tab:green"}[r['algorithm']]
        axs[0,idx].plot(r['power'],label=r['algorithm'],color=color)
        if r['algorithm']=="P&O":axs[0,idx].plot(r['mpp'],'k--',label="Theoretical MPP")
        axs[1,idx].plot(r['vpv'],label=r['algorithm'],color=color)
        summary['results'].append({"scenario":r['name'],"algorithm":r['algorithm'],"tracking_efficiency_pct":r['efficiency'],"final_power_w":float(r['power'][-1]),"final_mpp_w":float(r['mpp'][-1]),"final_vpv_v":float(r['vpv'][-1]),"final_vout_v":float(r['vout'][-1]),"last_duty":float(r['duty'][-1])})
    axs[0,0].set_title("Averaged boost: irradiance step");axs[0,1].set_title("Averaged boost: partial shading")
    for ax in axs[0,:]:ax.set_ylabel("PV input power (W)");ax.grid(True,alpha=.3);ax.legend(fontsize=8)
    for ax in axs[1,:]:ax.set_ylabel("PV voltage (V)");ax.set_xlabel("MPPT sample");ax.grid(True,alpha=.3);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/"averaged_boost_converter_comparison.png",dpi=180);plt.close(fig)
    summary['runtime_seconds']=time.time()-start
    (OUT/"averaged_boost_converter_results.json").write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print("Runtime seconds:",round(summary['runtime_seconds'],2))

if __name__=="__main__":main()
