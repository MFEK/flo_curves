use flo_curves::*;
use flo_curves::bezier::*;
use flo_curves::bezier::path::*;
use flo_curves::bezier::path::algorithms::*;

fn circle_ray_cast(circle_center: Coord2, radius: f64) -> impl Fn(Coord2, Coord2) -> Vec<RayCollision<Coord2, ()>> {
    move |from: Coord2, to: Coord2| {
        let from    = from - circle_center;
        let to      = to - circle_center;

        let x1      = from.x();
        let y1      = from.y();
        let x2      = to.x();
        let y2      = to.y();

        let dx      = x2-x1;
        let dy      = y2-y1;
        let dr      = (dx*dx + dy*dy).sqrt();

        let d       = x1*y2 - x2*y1;

        let xc1     = (d*dy + dy.signum()*dx)*((radius*radius*dr*dr - d*d).sqrt())/(dr*dr);
        let xc2     = (d*dy - dy.signum()*dx)*((radius*radius*dr*dr - d*d).sqrt())/(dr*dr);
        let yc1     = (-d*dx + dy.abs())*((radius*radius*dr*dr - d*d).sqrt())/(dr*dr);
        let yc2     = (-d*dx - dy.abs())*((radius*radius*dr*dr - d*d).sqrt())/(dr*dr);

        vec![
            RayCollision::new(Coord2(xc1, yc1)+circle_center, ()), RayCollision::new(Coord2(xc2, yc2)+circle_center, ())
        ]
    }
}

#[test]
fn ray_cast_to_circle() {
    let ray_cast    = circle_ray_cast(Coord2(10.0, 10.0), 5.0);

    let from_center = ray_cast(Coord2(10.0, 10.0), Coord2(11.0, 11.0));

    assert!(from_center.len() == 2);
    assert!((from_center[0].position.distance_to(&Coord2(10.0, 10.0))-5.0).abs() < 0.1);
    assert!((from_center[1].position.distance_to(&Coord2(10.0, 10.0))-5.0).abs() < 0.1);
    assert!(from_center[0].position.distance_to(&Coord2(13.54, 13.54)) < 0.1);
    assert!(from_center[1].position.distance_to(&Coord2(6.46, 6.46)) < 0.1);

    let offset = ray_cast(Coord2(11.0, 11.0), Coord2(12.0, 12.0));

    assert!(offset.len() == 2);
    assert!((offset[0].position.distance_to(&Coord2(10.0, 10.0))-5.0).abs() < 0.1);
    assert!((offset[1].position.distance_to(&Coord2(10.0, 10.0))-5.0).abs() < 0.1);
    assert!(offset[0].position.distance_to(&Coord2(13.54, 13.54)) < 0.1);
    assert!(offset[1].position.distance_to(&Coord2(6.46, 6.46)) < 0.1);

    let offset2 = ray_cast(Coord2(11.0, 11.0), Coord2(12.0, 11.0));

    assert!(offset2.len() == 2);
    assert!((offset2[0].position.distance_to(&Coord2(10.0, 10.0))-5.0).abs() < 0.1);
    assert!((offset2[1].position.distance_to(&Coord2(10.0, 10.0))-5.0).abs() < 0.1);
}

#[test]
fn fill_concave_circle() {
    // Simple circle ray-casting algorithm
    let circle_center   = Coord2(10.0, 10.0);
    let radius          = 5.0;
    let circle_ray_cast = circle_ray_cast(circle_center, radius);

    // Flood-fill this curve
    let path = flood_fill_concave::<SimpleBezierPath, _, _, _,_>(circle_center, &FillSettings::default(), circle_ray_cast);

    assert!(path.is_some());
    assert!(path.as_ref().unwrap().len() == 1);

    for curve in path.unwrap()[0].to_curves::<Curve<Coord2>>() {
        for t in 0..100 {
            let t           = (t as f64)/100.0;
            let distance    = circle_center.distance_to(&curve.point_at_pos(t));

            assert!((distance-radius).abs() < 1.0);
        }
    }
}