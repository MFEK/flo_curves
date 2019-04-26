use flo_curves::*;
use flo_curves::debug::*;
use flo_curves::bezier::path::*;

use super::svg::*;

#[test]
fn remove_interior_points_1() {
    // Complicated curve found in FlowBetween that produces 0 points when interior points are removed
    // It appears this has three curves that converge on a single point, which generates two points in the output, 
    // which in turn produces a spurious edge, which prevents us from being able to follow the path all the way around.
    let curve = BezierPathBuilder::<SimpleBezierPath>::start(Coord2(562.0692138671875, 669.944580078125))
        .curve_to((Coord2(562.0692138671875, 669.944580078125), Coord2(562.0692138671875, 669.944580078125)), Coord2(562.0692138671875, 669.944580078125))
        .curve_to((Coord2(562.4200439453125, 669.9562377929688), Coord2(562.6718139648438, 670.0160522460938)), Coord2(562.8291015625, 670.0160522460938))
        .curve_to((Coord2(562.7747802734375, 670.2525634765625), Coord2(563.1968383789063, 669.9431762695313)), Coord2(563.401611328125, 669.6962890625))
        .curve_to((Coord2(563.218505859375, 669.447021484375), Coord2(562.7468872070313, 668.9757690429688)), Coord2(562.6525268554688, 669.1633911132813))
        .curve_to((Coord2(562.3690185546875, 669.1181640625), Coord2(562.02490234375, 668.9761352539063)), Coord2(561.6610107421875, 668.9097290039063))
        .curve_to((Coord2(560.6327514648438, 668.3336181640625), Coord2(560.5078125, 668.680419921875)), Coord2(560.7962036132813, 668.913818359375))
        .curve_to((Coord2(560.7932739257813, 669.053955078125), Coord2(560.795166015625, 668.9855346679688)), Coord2(560.7711791992188, 669.1161499023438))
        .curve_to((Coord2(560.0169067382813, 670.2783813476563), Coord2(561.4442749023438, 669.9208984375)), Coord2(561.7700805664063, 668.3951416015625))
        .curve_to((Coord2(560.5978393554688, 668.2579956054688), Coord2(555.1843872070313, 665.8880615234375)), Coord2(552.6854248046875, 664.9908447265625))
        .curve_to((Coord2(552.62158203125, 664.9489135742188), Coord2(552.4188842773438, 664.8121948242188)), Coord2(552.1951293945313, 664.701171875))
        .curve_to((Coord2(552.0418090820313, 669.23193359375), Coord2(555.0795288085938, 667.7922973632813)), Coord2(555.9325561523438, 664.439697265625))
        .curve_to((Coord2(555.8035278320313, 663.3936767578125), Coord2(543.4547729492188, 664.566162109375)), Coord2(541.8832397460938, 667.4561767578125))
        .curve_to((Coord2(542.26611328125, 672.0006103515625), Coord2(548.4946899414063, 670.5872192382813)), Coord2(547.83984375, 666.5941772460938))
        .curve_to((Coord2(546.2003784179688, 665.4840087890625), Coord2(543.0369262695313, 665.3294677734375)), Coord2(543.1106567382813, 665.9275512695313))
        .curve_to((Coord2(536.3306274414063, 669.7837524414063), Coord2(541.8121337890625, 670.9800415039063)), Coord2(539.4649658203125, 666.6785888671875))
        .curve_to((Coord2(536.6891479492188, 665.93017578125), Coord2(534.8207397460938, 663.9938354492188)), Coord2(533.337890625, 661.9244995117188))
        .curve_to((Coord2(532.1223754882813, 662.1298828125), Coord2(530.9287109375, 662.1915893554688)), Coord2(534.033203125, 663.5484619140625))
        .curve_to((Coord2(539.8789672851563, 669.0048828125), Coord2(535.4338989257813, 664.3715209960938)), Coord2(530.1646118164063, 657.32666015625))
        .curve_to((Coord2(525.6614379882813, 654.2191162109375), Coord2(526.3388671875, 656.8445434570313)), Coord2(530.332275390625, 658.5115356445313))
        .curve_to((Coord2(530.9607543945313, 663.3235473632813), Coord2(535.1883544921875, 667.216552734375)), Coord2(533.1292724609375, 661.65673828125))
        .curve_to((Coord2(526.8078002929688, 654.7847290039063), Coord2(527.2481689453125, 655.82421875)), Coord2(528.5620727539063, 658.5321044921875))
        .curve_to((Coord2(529.048828125, 663.075927734375), Coord2(530.8765869140625, 662.1258544921875)), Coord2(531.5584106445313, 659.6661987304688))
        .curve_to((Coord2(530.1249389648438, 657.940185546875), Coord2(529.1561889648438, 657.2536010742188)), Coord2(528.7389526367188, 655.5059814453125))
        .curve_to((Coord2(527.8021240234375, 654.7122192382813), Coord2(529.899658203125, 656.5814819335938)), Coord2(531.8333740234375, 654.7963256835938))
        .curve_to((Coord2(538.0204467773438, 653.547119140625), Coord2(542.1532592773438, 652.2764892578125)), Coord2(544.957275390625, 652.1034545898438))
        .curve_to((Coord2(545.7479858398438, 652.0574340820313), Coord2(546.3248291015625, 651.8165283203125)), Coord2(546.8508911132813, 651.8157958984375))
        .curve_to((Coord2(548.2747802734375, 652.2127685546875), Coord2(548.1990356445313, 651.2047119140625)), Coord2(547.912109375, 650.8655395507813))
        .curve_to((Coord2(547.7791748046875, 650.193359375), Coord2(549.1414184570313, 650.476806640625)), Coord2(548.0958251953125, 650.5689086914063))
        .curve_to((Coord2(548.0958251953125, 650.7786865234375), Coord2(548.0958251953125, 651.0584716796875)), Coord2(548.0958251953125, 651.2682495117188))
        .curve_to((Coord2(548.9656982421875, 651.643798828125), Coord2(547.8914184570313, 651.3145141601563)), Coord2(549.1207275390625, 650.8655395507813))
        .curve_to((Coord2(548.8338012695313, 650.4108276367188), Coord2(547.700927734375, 649.2344360351563)), Coord2(546.8508911132813, 649.63134765625))
        .curve_to((Coord2(546.3272705078125, 649.630615234375), Coord2(545.5951538085938, 649.4255981445313)), Coord2(544.8019409179688, 649.4730834960938))
        .curve_to((Coord2(542.03857421875, 649.6287841796875), Coord2(537.2066040039063, 649.3989868164063)), Coord2(530.6641845703125, 650.7567138671875))
        .curve_to((Coord2(529.2568359375, 650.2000122070313), Coord2(525.3572998046875, 653.1232299804688)), Coord2(525.4530639648438, 656.3499145507813))
        .curve_to((Coord2(526.0859985351563, 658.6912231445313), Coord2(527.9020385742188, 660.7957763671875)), Coord2(529.1198120117188, 661.9281616210938))
        .curve_to((Coord2(532.0623779296875, 661.90576171875), Coord2(533.4664306640625, 660.0416259765625)), Coord2(529.3554077148438, 657.97998046875))
        .curve_to((Coord2(526.156005859375, 654.2037353515625), Coord2(522.1826782226563, 656.4036254882813)), Coord2(530.3896484375, 664.3438720703125))
        .curve_to((Coord2(536.754150390625, 667.3721923828125), Coord2(535.8456420898438, 660.3375854492188)), Coord2(530.91162109375, 658.0571899414063))
        .curve_to((Coord2(529.3756103515625, 652.6741333007813), Coord2(525.1596069335938, 652.45458984375)), Coord2(527.2052612304688, 659.278076171875))
        .curve_to((Coord2(532.2788696289063, 667.9177856445313), Coord2(540.3832397460938, 669.9564208984375)), Coord2(534.7332763671875, 662.905029296875))
        .curve_to((Coord2(532.432373046875, 658.3805541992188), Coord2(530.3565063476563, 660.47900390625)), Coord2(530.7684326171875, 663.4412841796875))
        .curve_to((Coord2(531.7405395507813, 665.5307006835938), Coord2(535.3882446289063, 669.1942138671875)), Coord2(538.6748046875, 669.9224853515625))
        .curve_to((Coord2(545.757080078125, 667.9179077148438), Coord2(541.6903686523438, 667.3967895507813)), Coord2(543.7351684570313, 669.8466186523438))
        .curve_to((Coord2(545.7384643554688, 670.13720703125), Coord2(545.3059692382813, 669.0546875)), Coord2(544.8681640625, 669.27587890625))
        .curve_to((Coord2(544.5274047851563, 665.6309204101563), Coord2(545.35498046875, 665.0091552734375)), Coord2(545.9674682617188, 668.1023559570313))
        .curve_to((Coord2(544.9426879882813, 667.5361328125), Coord2(553.4862670898438, 669.1529541015625)), Coord2(556.2109375, 667.8751831054688))
        .curve_to((Coord2(557.3668823242188, 664.4981079101563), Coord2(550.9618530273438, 662.6256103515625)), Coord2(549.96337890625, 666.1478271484375))
        .curve_to((Coord2(551.0449829101563, 667.5820922851563), Coord2(551.2767333984375, 667.6608276367188)), Coord2(551.4939575195313, 667.7685546875))
        .curve_to((Coord2(554.0316772460938, 668.719970703125), Coord2(560.2760620117188, 670.0292358398438)), Coord2(561.5628051757813, 670.177978515625))
        .curve_to((Coord2(562.9513549804688, 668.7757568359375), Coord2(560.3701782226563, 666.8861694335938)), Coord2(559.4100952148438, 668.6519775390625))
        .curve_to((Coord2(559.3812255859375, 668.7340698242188), Coord2(559.3759765625, 668.8223876953125)), Coord2(559.3749389648438, 668.913818359375))
        .curve_to((Coord2(559.663330078125, 669.5628662109375), Coord2(561.01806640625, 670.2659912109375)), Coord2(561.4695434570313, 669.9596557617188))
        .curve_to((Coord2(561.845458984375, 670.0281982421875), Coord2(562.2411499023438, 669.9994506835938)), Coord2(562.5125732421875, 670.0426025390625))
        .curve_to((Coord2(562.97314453125, 670.3184814453125), Coord2(562.8713989257813, 669.9105834960938)), Coord2(562.6882934570313, 669.6962890625))
        .curve_to((Coord2(562.89306640625, 669.4701538085938), Coord2(563.1842041015625, 669.1715698242188)), Coord2(562.8291015625, 669.4080810546875))
        .curve_to((Coord2(562.6779174804688, 669.4080810546875), Coord2(562.442626953125, 669.45654296875)), Coord2(562.085693359375, 669.44482421875))
        .build();

    // Create the graph path from the source side
    let mut merged_path = GraphPath::new();
    merged_path         = merged_path.merge(GraphPath::from_merged_paths(vec![&curve].into_iter().map(|path| (path, PathLabel(PathSource::Path1, PathDirection::from(path))))));

    // Collide the path with itself to find the intersections
    merged_path.self_collide(0.01);
    merged_path.set_exterior_by_removing_interior_points();
    merged_path.heal_exterior_gaps();
    println!("{}", graph_path_svg_string(&merged_path, vec![]));

    println!("{:?}", svg_path_string(&curve));
    let with_points_removed: Vec<SimpleBezierPath> = path_remove_interior_points(&vec![curve], 0.01);

    println!("{:?}", with_points_removed.iter()
        .map(|path| svg_path_string(path))
        .collect::<Vec<_>>());

    assert!(with_points_removed.len() > 0);
}

#[test]
fn remove_interior_points_1_without_failing_section() {
    // Complicated curve found in FlowBetween that produces 0 points when interior points are removed, variant with the section that was causing a failure removed
    let curve = BezierPathBuilder::<SimpleBezierPath>::start(Coord2(562.0692138671875, 669.944580078125))
        .curve_to((Coord2(562.0692138671875, 669.944580078125), Coord2(562.0692138671875, 669.944580078125)), Coord2(562.0692138671875, 669.944580078125))
        .curve_to((Coord2(562.4200439453125, 669.9562377929688), Coord2(562.6718139648438, 670.0160522460938)), Coord2(562.8291015625, 670.0160522460938))
        .curve_to((Coord2(562.7747802734375, 670.2525634765625), Coord2(563.1968383789063, 669.9431762695313)), Coord2(563.401611328125, 669.6962890625))
        .curve_to((Coord2(563.218505859375, 669.447021484375), Coord2(562.7468872070313, 668.9757690429688)), Coord2(562.6525268554688, 669.1633911132813))
        .curve_to((Coord2(562.3690185546875, 669.1181640625), Coord2(562.02490234375, 668.9761352539063)), Coord2(561.6610107421875, 668.9097290039063))
        .curve_to((Coord2(560.6327514648438, 668.3336181640625), Coord2(560.5078125, 668.680419921875)), Coord2(560.7962036132813, 668.913818359375))
        .curve_to((Coord2(560.7932739257813, 669.053955078125), Coord2(560.795166015625, 668.9855346679688)), Coord2(560.7711791992188, 669.1161499023438))
        .curve_to((Coord2(560.0169067382813, 670.2783813476563), Coord2(561.4442749023438, 669.9208984375)), Coord2(561.7700805664063, 668.3951416015625))
        .curve_to((Coord2(560.5978393554688, 668.2579956054688), Coord2(555.1843872070313, 665.8880615234375)), Coord2(552.6854248046875, 664.9908447265625))
        .curve_to((Coord2(552.62158203125, 664.9489135742188), Coord2(552.4188842773438, 664.8121948242188)), Coord2(552.1951293945313, 664.701171875))
        .curve_to((Coord2(552.0418090820313, 669.23193359375), Coord2(555.0795288085938, 667.7922973632813)), Coord2(555.9325561523438, 664.439697265625))
        .curve_to((Coord2(555.8035278320313, 663.3936767578125), Coord2(543.4547729492188, 664.566162109375)), Coord2(541.8832397460938, 667.4561767578125))
        .curve_to((Coord2(542.26611328125, 672.0006103515625), Coord2(548.4946899414063, 670.5872192382813)), Coord2(547.83984375, 666.5941772460938))
        .curve_to((Coord2(546.2003784179688, 665.4840087890625), Coord2(543.0369262695313, 665.3294677734375)), Coord2(543.1106567382813, 665.9275512695313))
        .curve_to((Coord2(536.3306274414063, 669.7837524414063), Coord2(541.8121337890625, 670.9800415039063)), Coord2(539.4649658203125, 666.6785888671875))
        .curve_to((Coord2(536.6891479492188, 665.93017578125), Coord2(534.8207397460938, 663.9938354492188)), Coord2(533.337890625, 661.9244995117188))
        .curve_to((Coord2(532.1223754882813, 662.1298828125), Coord2(530.9287109375, 662.1915893554688)), Coord2(534.033203125, 663.5484619140625))
        .curve_to((Coord2(539.8789672851563, 669.0048828125), Coord2(535.4338989257813, 664.3715209960938)), Coord2(530.1646118164063, 657.32666015625))
        .curve_to((Coord2(525.6614379882813, 654.2191162109375), Coord2(526.3388671875, 656.8445434570313)), Coord2(530.332275390625, 658.5115356445313))
        .curve_to((Coord2(530.9607543945313, 663.3235473632813), Coord2(535.1883544921875, 667.216552734375)), Coord2(533.1292724609375, 661.65673828125))
        .curve_to((Coord2(526.8078002929688, 654.7847290039063), Coord2(527.2481689453125, 655.82421875)), Coord2(528.5620727539063, 658.5321044921875))
        .curve_to((Coord2(529.048828125, 663.075927734375), Coord2(530.8765869140625, 662.1258544921875)), Coord2(531.5584106445313, 659.6661987304688))
        .curve_to((Coord2(530.1249389648438, 657.940185546875), Coord2(529.1561889648438, 657.2536010742188)), Coord2(528.7389526367188, 655.5059814453125))
        .curve_to((Coord2(527.8021240234375, 654.7122192382813), Coord2(529.899658203125, 656.5814819335938)), Coord2(531.8333740234375, 654.7963256835938))
        .curve_to((Coord2(538.0204467773438, 653.547119140625), Coord2(542.1532592773438, 652.2764892578125)), Coord2(544.957275390625, 652.1034545898438))
        .curve_to((Coord2(545.7479858398438, 652.0574340820313), Coord2(546.3248291015625, 651.8165283203125)), Coord2(546.8508911132813, 651.8157958984375))
        .curve_to((Coord2(548.2747802734375, 652.2127685546875), Coord2(548.1990356445313, 651.2047119140625)), Coord2(547.912109375, 650.8655395507813))
        .curve_to((Coord2(547.7791748046875, 650.193359375), Coord2(549.1414184570313, 650.476806640625)), Coord2(548.0958251953125, 650.5689086914063))
        .curve_to((Coord2(548.0958251953125, 650.7786865234375), Coord2(548.0958251953125, 651.0584716796875)), Coord2(548.0958251953125, 651.2682495117188))
        .curve_to((Coord2(548.9656982421875, 651.643798828125), Coord2(547.8914184570313, 651.3145141601563)), Coord2(549.1207275390625, 650.8655395507813))
        .curve_to((Coord2(548.8338012695313, 650.4108276367188), Coord2(547.700927734375, 649.2344360351563)), Coord2(546.8508911132813, 649.63134765625))
        .curve_to((Coord2(546.3272705078125, 649.630615234375), Coord2(545.5951538085938, 649.4255981445313)), Coord2(544.8019409179688, 649.4730834960938))
        .curve_to((Coord2(542.03857421875, 649.6287841796875), Coord2(537.2066040039063, 649.3989868164063)), Coord2(530.6641845703125, 650.7567138671875))
        .curve_to((Coord2(529.2568359375, 650.2000122070313), Coord2(525.3572998046875, 653.1232299804688)), Coord2(525.4530639648438, 656.3499145507813))
        .curve_to((Coord2(526.0859985351563, 658.6912231445313), Coord2(527.9020385742188, 660.7957763671875)), Coord2(529.1198120117188, 661.9281616210938))
        .curve_to((Coord2(532.0623779296875, 661.90576171875), Coord2(533.4664306640625, 660.0416259765625)), Coord2(529.3554077148438, 657.97998046875))
        .curve_to((Coord2(526.156005859375, 654.2037353515625), Coord2(522.1826782226563, 656.4036254882813)), Coord2(530.3896484375, 664.3438720703125))
        .curve_to((Coord2(536.754150390625, 667.3721923828125), Coord2(535.8456420898438, 660.3375854492188)), Coord2(530.91162109375, 658.0571899414063))
        .curve_to((Coord2(529.3756103515625, 652.6741333007813), Coord2(525.1596069335938, 652.45458984375)), Coord2(527.2052612304688, 659.278076171875))
        .curve_to((Coord2(532.2788696289063, 667.9177856445313), Coord2(540.3832397460938, 669.9564208984375)), Coord2(534.7332763671875, 662.905029296875))
        .curve_to((Coord2(532.432373046875, 658.3805541992188), Coord2(530.3565063476563, 660.47900390625)), Coord2(530.7684326171875, 663.4412841796875))
        /*FAIL*//*.curve_to((Coord2(531.7405395507813, 665.5307006835938), Coord2(535.3882446289063, 669.1942138671875)), Coord2(538.6748046875, 669.9224853515625))*/
        .curve_to((Coord2(545.757080078125, 667.9179077148438), Coord2(541.6903686523438, 667.3967895507813)), Coord2(543.7351684570313, 669.8466186523438))
        .curve_to((Coord2(545.7384643554688, 670.13720703125), Coord2(545.3059692382813, 669.0546875)), Coord2(544.8681640625, 669.27587890625))
        .curve_to((Coord2(544.5274047851563, 665.6309204101563), Coord2(545.35498046875, 665.0091552734375)), Coord2(545.9674682617188, 668.1023559570313))
        .curve_to((Coord2(544.9426879882813, 667.5361328125), Coord2(553.4862670898438, 669.1529541015625)), Coord2(556.2109375, 667.8751831054688))
        .curve_to((Coord2(557.3668823242188, 664.4981079101563), Coord2(550.9618530273438, 662.6256103515625)), Coord2(549.96337890625, 666.1478271484375))
        .curve_to((Coord2(551.0449829101563, 667.5820922851563), Coord2(551.2767333984375, 667.6608276367188)), Coord2(551.4939575195313, 667.7685546875))
        .curve_to((Coord2(554.0316772460938, 668.719970703125), Coord2(560.2760620117188, 670.0292358398438)), Coord2(561.5628051757813, 670.177978515625))
        .curve_to((Coord2(562.9513549804688, 668.7757568359375), Coord2(560.3701782226563, 666.8861694335938)), Coord2(559.4100952148438, 668.6519775390625))
        .curve_to((Coord2(559.3812255859375, 668.7340698242188), Coord2(559.3759765625, 668.8223876953125)), Coord2(559.3749389648438, 668.913818359375))
        .curve_to((Coord2(559.663330078125, 669.5628662109375), Coord2(561.01806640625, 670.2659912109375)), Coord2(561.4695434570313, 669.9596557617188))
        .curve_to((Coord2(561.845458984375, 670.0281982421875), Coord2(562.2411499023438, 669.9994506835938)), Coord2(562.5125732421875, 670.0426025390625))
        .curve_to((Coord2(562.97314453125, 670.3184814453125), Coord2(562.8713989257813, 669.9105834960938)), Coord2(562.6882934570313, 669.6962890625))
        .curve_to((Coord2(562.89306640625, 669.4701538085938), Coord2(563.1842041015625, 669.1715698242188)), Coord2(562.8291015625, 669.4080810546875))
        .curve_to((Coord2(562.6779174804688, 669.4080810546875), Coord2(562.442626953125, 669.45654296875)), Coord2(562.085693359375, 669.44482421875))
        .build();

    println!("{:?}", svg_path_string(&curve));
    let with_points_removed: Vec<SimpleBezierPath> = path_remove_interior_points(&vec![curve], 0.01);

    println!("{:?}", with_points_removed.iter()
        .map(|path| svg_path_string(path))
        .collect::<Vec<_>>());

    assert!(with_points_removed.len() > 0);
}

#[test]
fn remove_interior_points_2() {
    // Complicated curve found in FlowBetween that produces 0 points when interior points are removed
    // It appears this has three curves that converge on a single point, which generates two points in the output, 
    // which in turn produces a spurious edge, which prevents us from being able to follow the path all the way around.
    let curve = BezierPathBuilder::<SimpleBezierPath>::start(Coord2(589.8298950195313, 841.699951171875))
        .curve_to((Coord2(589.8298950195313, 841.699951171875), Coord2(589.8298950195313, 841.699951171875)), Coord2(589.8298950195313, 841.699951171875))
        .curve_to((Coord2(585.0781860351563, 841.545166015625), Coord2(588.116943359375, 846.1569213867188)), Coord2(589.9508056640625, 846.92041015625))
        .curve_to((Coord2(593.9074096679688, 850.3338623046875), Coord2(596.3680419921875, 855.8639526367188)), Coord2(600.2550048828125, 860.024169921875))
        .curve_to((Coord2(602.3019409179688, 864.72900390625), Coord2(603.487060546875, 861.721435546875)), Coord2(602.1428833007813, 859.0895385742188))
        .curve_to((Coord2(607.4638061523438, 858.4710693359375), Coord2(614.4444580078125, 855.14404296875)), Coord2(608.3931884765625, 855.6187133789063))
        .curve_to((Coord2(604.7843627929688, 851.9526977539063), Coord2(601.4735107421875, 847.9655151367188)), Coord2(597.78515625, 843.8760986328125))
        .curve_to((Coord2(601.0536499023438, 837.7391357421875), Coord2(590.90966796875, 841.439453125)), Coord2(587.8450927734375, 847.3414916992188))
        .curve_to((Coord2(592.2240600585938, 850.6311645507813), Coord2(595.8001098632813, 856.1324462890625)), Coord2(599.6971435546875, 861.4691772460938))
        .curve_to((Coord2(599.6600952148438, 866.2685546875), Coord2(601.5029907226563, 861.010498046875)), Coord2(601.408447265625, 857.6356811523438))
        .curve_to((Coord2(605.051025390625, 858.197509765625), Coord2(608.0866088867188, 854.1636352539063)), Coord2(597.3378295898438, 846.8604125976563))
        .curve_to((Coord2(597.2238159179688, 836.9576416015625), Coord2(590.7571411132813, 843.5430297851563)), Coord2(587.1199340820313, 848.599365234375))
        .curve_to((Coord2(588.7532348632813, 853.0540161132813), Coord2(591.633544921875, 856.119873046875)), Coord2(594.626708984375, 853.6188354492188))
        .curve_to((Coord2(596.7156982421875, 852.8362426757813), Coord2(595.0059814453125, 845.878662109375)), Coord2(591.52490234375, 845.5113525390625))
        .curve_to((Coord2(585.76171875, 847.6647338867188), Coord2(580.7750244140625, 855.853759765625)), Coord2(586.7627563476563, 853.3876342773438))
        .curve_to((Coord2(588.5208129882813, 859.3195190429688), Coord2(594.2566528320313, 860.6160278320313)), Coord2(592.3621826171875, 860.9254760742188))
        .curve_to((Coord2(594.9733276367188, 864.4375), Coord2(593.3421020507813, 848.7232055664063)), Coord2(586.76220703125, 847.8418579101563))
        .curve_to((Coord2(589.7845458984375, 841.6835327148438), Coord2(583.6079711914063, 848.498046875)), Coord2(580.9037475585938, 853.9146118164063))
        .curve_to((Coord2(580.701904296875, 853.186767578125), Coord2(578.50439453125, 857.2315063476563)), Coord2(581.5901489257813, 860.4940795898438))
        .curve_to((Coord2(585.6346435546875, 863.285400390625), Coord2(589.900146484375, 854.3807373046875)), Coord2(584.1525268554688, 856.2511596679688))
        .curve_to((Coord2(590.3831787109375, 852.05712890625), Coord2(578.9157104492188, 850.2012329101563)), Coord2(574.5430297851563, 856.5203247070313))
        .curve_to((Coord2(573.6943969726563, 863.1355590820313), Coord2(580.0052490234375, 871.26220703125)), Coord2(575.3004760742188, 871.1060791015625))
        .curve_to((Coord2(576.81103515625, 870.624267578125), Coord2(572.30712890625, 859.2913818359375)), Coord2(570.9198608398438, 861.718994140625))
        .curve_to((Coord2(572.5287475585938, 864.7382202148438), Coord2(581.41259765625, 882.9050903320313)), Coord2(580.4722900390625, 881.7498779296875))
        .curve_to((Coord2(580.0606689453125, 880.2344970703125), Coord2(575.6553955078125, 869.0311889648438)), Coord2(573.716552734375, 868.6065673828125))
        .curve_to((Coord2(570.4192504882813, 866.5391845703125), Coord2(572.1432495117188, 889.7837524414063)), Coord2(575.9349365234375, 889.2540893554688))
        .curve_to((Coord2(579.9112548828125, 889.1182250976563), Coord2(573.3362426757813, 870.1537475585938)), Coord2(570.325439453125, 872.933349609375))
        .curve_to((Coord2(566.7039184570313, 872.4866333007813), Coord2(575.889892578125, 896.3516845703125)), Coord2(580.193359375, 885.1004028320313))
        .curve_to((Coord2(578.9361572265625, 882.8379516601563), Coord2(578.29638671875, 880.9623413085938)), Coord2(577.2049560546875, 878.0570678710938))
        .curve_to((Coord2(576.3244018554688, 875.5227661132813), Coord2(575.8396606445313, 874.0106811523438)), Coord2(575.3523559570313, 871.5857543945313))
        .curve_to((Coord2(567.6146240234375, 879.8153076171875), Coord2(569.26904296875, 890.168212890625)), Coord2(572.8831176757813, 890.166259765625))
        .curve_to((Coord2(580.7759399414063, 887.835693359375), Coord2(580.0247802734375, 885.56103515625)), Coord2(572.6173095703125, 889.1515502929688))
        .curve_to((Coord2(572.6390991210938, 889.1546020507813), Coord2(572.2571411132813, 889.167724609375)), Coord2(572.2820434570313, 889.1630249023438))
        .curve_to((Coord2(570.7896728515625, 887.8728637695313), Coord2(567.4065551757813, 888.0462036132813)), Coord2(572.1813354492188, 892.5457763671875))
        .curve_to((Coord2(570.9942016601563, 894.4725341796875), Coord2(577.9598999023438, 900.7188720703125)), Coord2(582.4383544921875, 902.0015258789063))
        .curve_to((Coord2(582.8182373046875, 902.308349609375), Coord2(586.3283081054688, 901.2371826171875)), Coord2(586.35205078125, 900.798583984375))
        .curve_to((Coord2(588.947998046875, 898.2053833007813), Coord2(592.195068359375, 891.016845703125)), Coord2(591.5047607421875, 889.0786743164063))
        .curve_to((Coord2(592.836669921875, 884.303955078125), Coord2(592.759033203125, 882.3919677734375)), Coord2(593.544921875, 881.51806640625))
        .curve_to((Coord2(594.1064453125, 880.7155151367188), Coord2(593.8582153320313, 881.4864501953125)), Coord2(596.4064331054688, 879.8722534179688))
        .curve_to((Coord2(597.3624877929688, 879.4691162109375), Coord2(597.849365234375, 879.2901611328125)), Coord2(598.5863037109375, 879.035400390625))
        .curve_to((Coord2(598.9070434570313, 878.928466796875), Coord2(599.0929565429688, 878.8623657226563)), Coord2(599.3098754882813, 878.7935180664063))
        .curve_to((Coord2(596.8707275390625, 882.4271240234375), Coord2(601.0760498046875, 876.7950439453125)), Coord2(603.7096557617188, 873.0940551757813))
        .curve_to((Coord2(603.7099609375, 873.09423828125), Coord2(603.7090454101563, 872.9913940429688)), Coord2(603.708740234375, 872.9912109375))
        .curve_to((Coord2(602.2408447265625, 869.050048828125), Coord2(594.5162353515625, 860.6947021484375)), Coord2(590.5521850585938, 859.4874267578125))
        .curve_to((Coord2(584.2527465820313, 852.785888671875), Coord2(581.061279296875, 852.8796997070313)), Coord2(581.9452514648438, 853.1057739257813))
        .curve_to((Coord2(582.593017578125, 852.9660034179688), Coord2(580.7717895507813, 856.9833984375)), Coord2(580.3909301757813, 856.538818359375))
        .curve_to((Coord2(580.433837890625, 856.8338623046875), Coord2(577.17236328125, 858.1321411132813)), Coord2(578.2269897460938, 857.4826049804688))
        .curve_to((Coord2(579.6019897460938, 857.3278198242188), Coord2(581.242431640625, 858.7636108398438)), Coord2(586.4614868164063, 860.9908447265625))
        .curve_to((Coord2(589.0213012695313, 862.7692260742188), Coord2(595.4768676757813, 865.4718017578125)), Coord2(598.6329345703125, 865.7658081054688))
        .curve_to((Coord2(598.567138671875, 866.820556640625), Coord2(603.420166015625, 864.4375610351563)), Coord2(603.9707641601563, 863.19189453125))
        .curve_to((Coord2(604.5771484375, 862.9888916015625), Coord2(605.8325805664063, 859.209716796875)), Coord2(605.4430541992188, 858.7909545898438))
        .curve_to((Coord2(604.993408203125, 855.49951171875), Coord2(601.7562866210938, 849.7847900390625)), Coord2(600.6087036132813, 850.0838623046875))
        .curve_to((Coord2(598.4024047851563, 846.3101196289063), Coord2(598.4458618164063, 848.7549438476563)), Coord2(599.054931640625, 849.7504272460938))
        .curve_to((Coord2(599.3753051757813, 849.5570678710938), Coord2(598.5122680664063, 852.48095703125)), Coord2(598.3043212890625, 852.2940063476563))
        .curve_to((Coord2(594.026123046875, 853.1077270507813), Coord2(590.7466430664063, 857.636474609375)), Coord2(597.5106811523438, 857.97314453125))
        .curve_to((Coord2(599.0369262695313, 859.3222045898438), Coord2(599.9699096679688, 856.8018798828125)), Coord2(600.7046508789063, 857.8336181640625))
        .curve_to((Coord2(602.2219848632813, 857.0648803710938), Coord2(604.2450561523438, 856.0399780273438)), Coord2(605.7623901367188, 855.271240234375))
        .curve_to((Coord2(606.0952758789063, 855.509765625), Coord2(607.6958618164063, 852.7109985351563)), Coord2(605.1021118164063, 850.0916748046875))
        .curve_to((Coord2(607.8820190429688, 846.5908203125), Coord2(593.875244140625, 843.7130737304688)), Coord2(588.6094360351563, 846.0635986328125))
        .curve_to((Coord2(588.1766357421875, 846.2265625), Coord2(587.21044921875, 849.5330200195313)), Coord2(587.5308227539063, 849.7504272460938))
        .curve_to((Coord2(588.139892578125, 853.0325317382813), Coord2(591.3778076171875, 858.6402587890625)), Coord2(592.365966796875, 858.1364135742188))
        .curve_to((Coord2(594.4129028320313, 861.7054443359375), Coord2(594.3702392578125, 859.3676147460938)), Coord2(593.9205932617188, 858.7909545898438))
        .curve_to((Coord2(593.531005859375, 859.0397338867188), Coord2(594.5933837890625, 855.8880004882813)), Coord2(594.7660522460938, 856.260986328125))
        .curve_to((Coord2(595.0291137695313, 855.397216796875), Coord2(599.3807983398438, 853.1875)), Coord2(598.6329345703125, 854.2422485351563))
        .curve_to((Coord2(598.4532470703125, 854.5361938476563), Coord2(597.2518920898438, 853.0940551757813)), Coord2(591.7156982421875, 850.7276611328125))
        .curve_to((Coord2(588.8387451171875, 848.8101196289063), Coord2(581.9440307617188, 846.1011962890625)), Coord2(578.2269897460938, 845.9464111328125))
        .curve_to((Coord2(577.8850708007813, 845.296875), Coord2(573.4860229492188, 846.9068603515625)), Coord2(572.7304077148438, 847.910888671875))
        .curve_to((Coord2(571.8158569335938, 847.940185546875), Coord2(569.74658203125, 852.5507202148438)), Coord2(570.394287109375, 853.1057739257813))
        .curve_to((Coord2(571.2783203125, 857.921142578125), Coord2(578.909423828125, 867.0386962890625)), Coord2(583.4408569335938, 868.6987915039063))
        .curve_to((Coord2(590.3077392578125, 875.8531494140625), Coord2(593.4223022460938, 875.197265625)), Coord2(591.9874877929688, 873.194091796875))
        .curve_to((Coord2(591.9872436523438, 873.1967163085938), Coord2(591.986328125, 873.0939331054688)), Coord2(591.9866333007813, 873.0940551757813))
        .curve_to((Coord2(594.6202392578125, 869.4052734375), Coord2(598.6082153320313, 863.8417358398438)), Coord2(595.7800903320313, 867.597900390625))
        .curve_to((Coord2(595.6080322265625, 867.6517333984375), Coord2(595.2333374023438, 867.7623901367188)), Coord2(594.86767578125, 867.8843994140625))
        .curve_to((Coord2(594.2318115234375, 868.0874633789063), Coord2(592.8426513671875, 868.5746459960938)), Coord2(591.7855224609375, 869.0294189453125))
        .curve_to((Coord2(590.3074340820313, 869.1311645507813), Coord2(585.4846801757813, 872.3850708007813)), Coord2(583.8530883789063, 874.7000122070313))
        .curve_to((Coord2(582.0348510742188, 877.52783203125), Coord2(580.2197875976563, 883.7698364257813)), Coord2(579.9202880859375, 886.5905151367188))
        .curve_to((Coord2(577.5986328125, 892.2477416992188), Coord2(579.3203735351563, 892.0960083007813)), Coord2(579.736328125, 890.97021484375))
        .curve_to((Coord2(579.3685302734375, 890.795166015625), Coord2(582.47412109375, 889.8472900390625)), Coord2(582.4383544921875, 890.154052734375))
        .curve_to((Coord2(583.9168090820313, 891.4367065429688), Coord2(587.2955322265625, 891.2642211914063)), Coord2(582.5258178710938, 886.7721557617188))
        .curve_to((Coord2(583.718017578125, 884.8529663085938), Coord2(576.7567138671875, 878.6074829101563)), Coord2(572.2820434570313, 877.3173217773438))
        .curve_to((Coord2(572.260009765625, 877.3126220703125), Coord2(571.8314208984375, 877.3272705078125)), Coord2(571.8067016601563, 877.3335571289063))
        .curve_to((Coord2(560.010498046875, 881.22509765625), Coord2(569.1023559570313, 904.441650390625)), Coord2(580.4130249023438, 899.290283203125))
        .curve_to((Coord2(587.358642578125, 896.5389404296875), Coord2(577.4598388671875, 865.2567138671875)), Coord2(569.1909790039063, 866.8143920898438))
        .curve_to((Coord2(564.2433471679688, 876.4852905273438), Coord2(565.2119750976563, 879.7418823242188)), Coord2(566.0805053710938, 882.063720703125))
        .curve_to((Coord2(566.91650390625, 884.5101928710938), Coord2(568.3486328125, 888.1217651367188)), Coord2(569.8612060546875, 890.8428344726563))
        .curve_to((Coord2(587.60546875, 903.775146484375), Coord2(589.3789672851563, 862.3728637695313)), Coord2(570.9526977539063, 861.138671875))
        .curve_to((Coord2(555.2870483398438, 863.2453002929688), Coord2(565.8951416015625, 904.8345336914063)), Coord2(580.59521484375, 900.0840454101563))
        .curve_to((Coord2(594.0852661132813, 895.3810424804688), Coord2(579.9393920898438, 852.4515380859375)), Coord2(565.9549560546875, 859.7515869140625))
        .curve_to((Coord2(557.6361083984375, 864.9191284179688), Coord2(571.5972900390625, 896.5897216796875)), Coord2(582.2666625976563, 893.362060546875))
        .curve_to((Coord2(596.2720947265625, 889.8972778320313), Coord2(586.6115112304688, 848.4343872070313)), Coord2(571.9376220703125, 850.0352172851563))
        .curve_to((Coord2(556.0201416015625, 851.1971435546875), Coord2(561.6185302734375, 885.7030029296875)), Coord2(577.4126586914063, 882.59521484375))
        .curve_to((Coord2(589.552978515625, 879.34228515625), Coord2(591.1780395507813, 853.8557739257813)), Coord2(584.5911254882813, 850.6495971679688))
        .curve_to((Coord2(573.9320678710938, 846.2090454101563), Coord2(564.6778564453125, 858.427490234375)), Coord2(573.91796875, 861.5853881835938))
        .curve_to((Coord2(572.7061767578125, 870.7952880859375), Coord2(589.696533203125, 873.0230712890625)), Coord2(593.0504760742188, 859.9937133789063))
        .curve_to((Coord2(595.9219970703125, 858.3511352539063), Coord2(586.99267578125, 843.552978515625)), Coord2(581.2372436523438, 842.6605834960938))
        .curve_to((Coord2(576.9470825195313, 848.0301513671875), Coord2(573.3963623046875, 858.2985229492188)), Coord2(577.3682861328125, 853.6705322265625))
        .curve_to((Coord2(573.3412475585938, 856.9037475585938), Coord2(593.6327514648438, 872.11083984375)), Coord2(601.6283569335938, 866.2468872070313))
        .curve_to((Coord2(603.8673706054688, 859.3587036132813), Coord2(597.4412231445313, 846.4760131835938)), Coord2(595.3387451171875, 847.1231689453125))
        .curve_to((Coord2(600.8814697265625, 844.0478515625), Coord2(586.74169921875, 842.0896606445313)), Coord2(580.974609375, 845.4783935546875))
        .curve_to((Coord2(577.481201171875, 849.10498046875), Coord2(595.04345703125, 864.8500366210938)), Coord2(599.9298095703125, 862.3732299804688))
        .curve_to((Coord2(603.3001708984375, 859.6436767578125), Coord2(598.7106323242188, 838.4157104492188)), Coord2(590.7146606445313, 839.1911010742188))
        .curve_to((Coord2(586.1773071289063, 843.9035034179688), Coord2(581.3427124023438, 853.1783447265625)), Coord2(589.6311645507813, 853.3121948242188))
        .curve_to((Coord2(591.4857788085938, 861.0636596679688), Coord2(603.1357421875, 865.7894287109375)), Coord2(608.3871459960938, 864.779541015625))
        .curve_to((Coord2(609.6311645507813, 860.09716796875), Coord2(608.9767456054688, 852.512939453125)), Coord2(607.8642578125, 855.7709350585938))
        .curve_to((Coord2(604.8263549804688, 851.1680908203125), Coord2(599.0707397460938, 843.4915161132813)), Coord2(593.83251953125, 839.4678955078125))
        .curve_to((Coord2(588.12841796875, 843.3626708984375), Coord2(584.9580688476563, 854.13720703125)), Coord2(590.4581909179688, 850.477294921875))
        .curve_to((Coord2(593.902099609375, 854.3041381835938), Coord2(597.5442504882813, 858.5637817382813)), Coord2(601.3976440429688, 862.492431640625))
        .curve_to((Coord2(596.497314453125, 864.138427734375), Coord2(605.6493530273438, 863.7943115234375)), Coord2(611.2799682617188, 862.292236328125))
        .curve_to((Coord2(610.607666015625, 857.74365234375), Coord2(606.7666625976563, 851.5074462890625)), Coord2(606.7361450195313, 853.9833984375))
        .curve_to((Coord2(602.8890380859375, 849.8455200195313), Coord2(597.0514526367188, 846.7515869140625)), Coord2(593.0549926757813, 843.3157958984375))
        .curve_to((Coord2(591.688232421875, 841.3230590820313), Coord2(585.3775024414063, 841.3017578125)), Coord2(589.6620483398438, 842.685791015625))
        .build();

    // Create the graph path from the source side
    let mut merged_path = GraphPath::new();
    merged_path         = merged_path.merge(GraphPath::from_merged_paths(vec![&curve].into_iter().map(|path| (path, PathLabel(PathSource::Path1, PathDirection::from(path))))));

    // Collide the path with itself to find the intersections
    merged_path.self_collide(0.01);
    merged_path.set_exterior_by_removing_interior_points();
    merged_path.heal_exterior_gaps();
    println!("{}", graph_path_svg_string(&merged_path, vec![]));

    println!("{:?}", svg_path_string(&curve));
    let with_points_removed: Vec<SimpleBezierPath> = path_remove_interior_points(&vec![curve], 0.01);

    println!("{:?}", with_points_removed.iter()
        .map(|path| svg_path_string(path))
        .collect::<Vec<_>>());

    assert!(with_points_removed.len() > 0);
}
